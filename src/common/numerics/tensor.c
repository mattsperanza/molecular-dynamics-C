// Author(s): Matthew Speranza
#include "../include/tensor.h"

#include <assert.h>
#include <math.h>
#include <string.h>

// Many of these methods are taken from FFX and are not optimized for performance.
// Performance is gained by using "code that writes code" so there is no looping/array overhead.

long factorial(const int n) {
  long val = 1;
  for(int i = 1; i < n; i++) {
    val += val*i;
  }
  return val;
}

int nChooseK(int n, int k) {
  return factorial(n) / (factorial(k) * factorial(n-k));
}

/**
 *
 * @param src source terms to be filled
 * @param r distance between i and j
 * @param beta ewald beta
 * @param tensorOrder
 */
void ewaldSource(REAL* src, REAL* r, REAL beta, int tensorOrder) {
  // Calculate complementary Boys function terms from Sagui et al. Eq. 2.24
  REAL prefactor = 2.0 * beta / sqrt(M_PI);
  REAL twoBetaSq = -2.0 * beta * beta;
  REAL R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  REAL betaR = beta * R;
  REAL betaRSq = betaR * betaR;
  REAL i2BetaRSq = 1 / (2.0*betaRSq);
  REAL expBRSq = exp(-betaRSq);
  double Fnc = sqrt(M_PI) * erfc(betaR) / (2.0 * betaR);
  for(int i = 0; i <= tensorOrder; i++) {
    src[i] = prefactor * pow(twoBetaSq, i);
    src[i] = src[i] * Fnc;
    Fnc = fma(fma(2.0, i, 1.0), Fnc, expBRSq) * i2BetaRSq;
  }
}

void tholeSource(REAL* src, REAL* r, REAL ai, REAL aj) {

}

/**
 * Something I actually did come up with.
 *
 * Better indexing for the tensor to access elements in
 * "tensor order" order so elements are not accessed essentially randomly.
 *
 * @param dx order x derivative
 * @param dy order y derivative
 * @param dz order z derivative
 * @param order tensor order
 */
int ti(int dx, int dy, int dz, int order) {
  int sumTo[order];
  sumTo[0] = 1;
  for(int i = 1; i < order; i++){
    sumTo[i] = sumTo[i-1] + i + 1;
  }
  int sumToIndex = dx + dy + dz;
  int index = 0;
  for(int i = 0; i < sumToIndex; i++) {
    index += sumTo[i];
  }
  for(int i = 0; i < dz; i++) {
    index += sumToIndex + 1;
    sumToIndex--;
  }
  index += dy;
  return index;
}

/**
 * Turns index into derivative orders.
 * @param index index to reverse
 * @param order
 * @return pointer to array of derivative orders dx, dy, dz
 */
void it(const int index, const int order, int* ret) {
  int sumTo[order+1];
  sumTo[0] = 1;
  for(int i = 1; i < order+1; i++){
    sumTo[i] = sumTo[i-1] + i + 1;
  }
  int val = index;
  int count = 0;
  while(val >= 0) {
    val -= sumTo[count];
    count++;
  }
  int sumToIndex = count-1;
  val = index;
  for(int i = 0; i < sumToIndex; i++) {
    val -= sumTo[i];
  }
  int sumToTemp = sumToIndex;
  count = 0;
  while(val >= 0) {
    val -= sumToTemp + 1;
    sumToTemp--;
    count++;
  }
  int dz = count-1;
  int dy = val + sumToTemp+2;
  int dx = sumToIndex - dy - dz;
  ret[0] = dx;
  ret[1] = dy;
  ret[2] = dz;
}

/*
Tlmnj recursive method from FFX - theory in paper more clear here than unrolled loop version:
if (m == 0 && n == 0) {
  if (l > 1) {
    return r[0] * Tlmnj(l - 1, 0, 0, j + 1, r, T000)
        + (l - 1) * Tlmnj(l - 2, 0, 0, j + 1, r, T000);
  } else if (l == 1) { // l == 1; d/dx is done.
    return r[0] * Tlmnj(0, 0, 0, j + 1, r, T000);
  } else {
    // l = m = n = 0; Recursion is done.
    return T000[j];
  }
} else if (n == 0) {
  // m >= 1
  if (m > 1) {
    return r[1] * Tlmnj(l, m - 1, 0, j + 1, r, T000)
        + (m - 1) * Tlmnj(l, m - 2, 0, j + 1, r, T000);
  }
  // m == 1; d/dy is done.
  return r[1] * Tlmnj(l, 0, 0, j + 1, r, T000);
} else {
  // n >= 1
  if (n > 1) {
    return r[2] * Tlmnj(l, m, n - 1, j + 1, r, T000)
        + (n - 1) * Tlmnj(l, m, n - 2, j + 1, r, T000);
  }
  // n == 1; d/dz is done.
  return r[2] * Tlmnj(l, m, 0, j + 1, r, T000);
}
*/

/**
 * Straight copy and paste from FFX CoulombTensorGlobal.java.
 *
 * @param tensor tensor to fill
 * @param r distance between i and j
 * @param src "source" terms for this tensor
 * @param order
 */
void generateTensor(REAL* tensor, REAL* r, REAL* srcTerm, int order) {
   const REAL x = r[0];
   const REAL y = r[1];
   const REAL z = r[2];
   const int o1 = order + 1;
   const int im = o1 * o1;
   const int in = im*o1;
   REAL storeTerms[in*o1];
   memcpy(storeTerms, srcTerm, (order+1)*sizeof(REAL));
   tensor[0] = storeTerms[0];
    // Find (d/dx)^l for l = 1..order (m = 0, n = 0)
    // Any (d/dx) term can be formed as
    // Tl00j = x * T(l-1)00(j+1) + (l-1) * T(l-2)00(j+1)
    // All intermediate terms are indexed as l*o1 + m*im + n*in + j;
    double current;
    double previous = storeTerms[1];
    // Store the l=1 tensor T100 (d/dx)
    tensor[ti(1, 0, 0, order)] = x * previous;
    // Starting the loop at l=2 avoids an if statement.
    for (int l = 2; l < o1; l++) {
      // Initial condition for the inner loop is formation of T100(l-1).
      // Starting the inner loop at a=1 avoids an if statement.
      // T100(l-1) = x * T000(l)
      current = x * storeTerms[l];
      int iw = o1 + l - 1;
      storeTerms[iw] = current;
      for (int a = 1; a < l - 1; a++) {
        // T200(l-2) = x * T100(l-1) + (2 - 1) * T000(l-1)
        // T300(l-3) = x * T200(l-2) + (3 - 1) * T100(l-2)
        // ...
        // T(l-1)001 = x * T(l-2)002 + (l - 2) * T(l-3)002
        current = fma(x, current, a * storeTerms[iw - o1]);
        iw += o1 - 1;
        storeTerms[iw] = current;
      }
      // Store the Tl00 tensor (d/dx)^l
      // Tl00 = x * T(l-1)001 + (l - 1) * T(l-2)001
      tensor[ti(l, 0, 0, order)] = fma(x, current, (l - 1) * previous);
      previous = current;
    }
    // Find (d/dx)^l * (d/dy)^m for l+m = 1..order (m > 0, n = 0)
    // Any (d/dy) term can be formed as:
    // Tlm0j = y * Tl(m-1)00(j+1) + (m-1) * Tl(m-2)00(j+1)
    for (int l = 0; l < order; l++) {
      // Store the m=1 tensor (d/dx)^l *(d/dy)
      // Tl10 = y * Tl001
      previous = storeTerms[l * o1 + 1];
      tensor[ti(l, 1, 0, order)] = y * previous;
      for (int m = 2; m + l < o1; m++) {
        // Tl10(m-1) = y * Tl00m;
        int iw = l * o1 + m;
        current = y * storeTerms[iw];
        iw += im - 1;
        storeTerms[iw] = current;
        for (int a = 1; a < m - 1; a++) {
          // Tl20(m-2) = Y * Tl10(m-1) + (2 - 1) * T100(m-1)
          // Tl30(m-3) = Y * Tl20(m-2) + (3 - 1) * Tl10(m-2)
          // ...
          // Tl(m-1)01 = Y * Tl(m-2)02 + (m - 2) * T(m-3)02
          current = fma(y, current, a * storeTerms[iw - im]);
          iw += im - 1;
          storeTerms[iw] = current;
        }
        // Store the tensor (d/dx)^l * (d/dy)^m
        // Tlm0 = y * Tl(m-1)01 + (m - 1) * Tl(m-2)01
        tensor[ti(l, m, 0, order)] = fma(y, current, (m - 1) * previous);
        previous = current;
      }
    }
    // Find (d/dx)^l * (d/dy)^m * (d/dz)^n for l+m+n = 1..order (n > 0)
    // Any (d/dz) term can be formed as:
    // Tlmnj = z * Tlm(n-1)(j+1) + (n-1) * Tlm(n-2)(j+1)
    for (int l = 0; l < order; l++) {
      for (int m = 0; m + l < order; m++) {
        // Store the n=1 tensor (d/dx)^l *(d/dy)^m * (d/dz)
        // Tlmn = z * Tlm01
        const int lm = m + l;
        const int lo1mim = l * o1 + m * im;
        previous = storeTerms[lo1mim + 1];
        tensor[ti(l, m, 1, order)] = z * previous;
        for (int n = 2; lm + n < o1; n++) {
          // Tlm1(n-1) = z * Tlm0n;
          int iw = lo1mim + n;
          current = z * storeTerms[iw];
          iw += in - 1;
          storeTerms[iw] = current;
          const int n1 = n - 1;
          for (int a = 1; a < n1; a++) {
            // Tlm2(n-2) = z * Tlm1(n-1) + (2 - 1) * T1m0(n-1)
            // Tlm3(n-3) = z * Tlm2(n-2) + (3 - 1) * Tlm1(n-2)
            // ...
            // Tlm(n-1)1 = z * Tlm(n-2)2 + (n - 2) * Tlm(n-3)2
            current = fma(z, current,  a * storeTerms[iw - in]);
            iw += in - 1;
            storeTerms[iw] = current;
          }
          // Store the tensor (d/dx)^l * (d/dy)^m * (d/dz)^n
          // Tlmn = z * Tlm(n-1)1 + (n - 1) * Tlm(n-2)1
          tensor[ti(l, m, n, order)] = fma(z, current, n1 * previous);
          previous = current;
        }
      }
    }
}

/**
 * Not for user use.
 * @param mpole
 * @param tensor
 * @param l
 * @param m
 * @param n
 * @param order
 * @return
 */
REAL contractMultipoleI(const REAL* mpole, const REAL* tensor, const int l, const int m, const int n, const int order) {
  REAL total = 0.0;
  // Charge
  total += mpole[0] * tensor[ti(l, m, n, order)];
  // Dipole
  total -= mpole[1] * tensor[ti(l+1, m, n, order)];
  total -= mpole[2] * tensor[ti(l, m+1, n, order)];
  total -= mpole[3] * tensor[ti(l, m, n+1, order)];
  // Trace
  total += mpole[4] * tensor[ti(l+2, m, n, order)];
  total += mpole[5] * tensor[ti(l, m+2, n, order)];
  total += mpole[6] * tensor[ti(l, m, n+2, order)];
  // Off-diagonal
  total += mpole[7] * tensor[ti(l+1, m+1, n, order)]; // xy
  total += mpole[8] * tensor[ti(l+1, m, n+1, order)]; // xz
  total += mpole[9] * tensor[ti(l, m+1, n+1, order)]; // yz
  return total;
}

/**
 * Not for user use.
 * @param mpole
 * @param tensor
 * @param l
 * @param m
 * @param n
 * @param order
 * @return
 */
REAL contractMultipoleJ(const REAL* mpole, const REAL* tensor, const int l, const int m, const int n, const int order) {
  REAL total = 0.0;
  // Charge
  total += mpole[0] * tensor[ti(l, m, n, order)];
  // Dipole
  total += mpole[1] * tensor[ti(l+1, m, n, order)]; // x
  total += mpole[2] * tensor[ti(l, m+1, n, order)]; // y
  total += mpole[3] * tensor[ti(l, m, n+1, order)]; // z
  // Trace
  total += mpole[4] * tensor[ti(l+2, m, n, order)]; // xx
  total += mpole[5] * tensor[ti(l, m+2, n, order)]; // yy
  total += mpole[6] * tensor[ti(l, m, n+2, order)]; // zz
  // Off-diagonal
  total += mpole[7] * tensor[ti(l+1, m+1, n, order)]; // xy
  total += mpole[8] * tensor[ti(l+1, m, n+1, order)]; // xz
  total += mpole[9] * tensor[ti(l, m+1, n+1, order)]; // yz
  return total;
}

/**
 * Fills E field components due to multipole I.
 * @param E electric field components
 * @param fieldOrder derivative order of electric field
 * @param mpole multipole
 * @param iOrJ true if I, false if J
 * @param tensor tensor
 * @param tensorOrder
 */
void potentialDueToMultipole(REAL* E, const int fieldOrder, const REAL* mpole, const bool iOrJ, const REAL* tensor, const int tensorOrder) {
  assert(fieldOrder <= tensorOrder-2);
  // There are fieldOrder+3 Choose 3 terms in the electric field
  int fieldSize = nChooseK(fieldOrder+3, 3);
  for(int i = 0; i < fieldSize; i++) {
    int derivatives[3] = {0, 0, 0};
    it(i, fieldOrder, derivatives);
    E[i] = iOrJ ? contractMultipoleI(mpole, tensor, derivatives[0], derivatives[1], derivatives[2], tensorOrder) :
    contractMultipoleJ(mpole, tensor, derivatives[0], derivatives[1], derivatives[2], tensorOrder);
  }
}

REAL multipoleEnergy(const REAL* mpole, const REAL* E) {
  // Indexing avoids calling ti()
  double total = mpole[0]*E[0];
  total = fma(mpole[1], E[1], total);
  total = fma(mpole[2], E[2], total);
  total = fma(mpole[3], E[3], total);
  total = fma(mpole[4], E[4], total);
  total = fma(mpole[5], E[6], total);
  total = fma(mpole[6], E[9], total);
  total = fma(mpole[7], E[5], total);
  total = fma(mpole[8], E[7], total);
  total = fma(mpole[9], E[8], total);
  return total;
}

REAL multipoleFieldInteraction(const REAL* mpole, const REAL* E, const int l, const int m, const int n, const int fieldOrder) {
  double total = mpole[0]*E[ti(l,m,n,fieldOrder)];
  // Dipole
  total = fma(mpole[1], E[ti(l+1,m,n,fieldOrder)], total);
  total = fma(mpole[2], E[ti(l,m+1,n,fieldOrder)], total);
  total = fma(mpole[3], E[ti(l,m,n+1,fieldOrder)], total);
  // Trace
  total = fma(mpole[4], E[ti(l+2,m,n,fieldOrder)], total);
  total = fma(mpole[5], E[ti(l,m+2,n,fieldOrder)], total);
  total = fma(mpole[6], E[ti(l,m,n+2,fieldOrder)], total);
  // Off-diagonal
  total = fma(mpole[7], E[ti(l+1,m+1,n,fieldOrder)], total);
  total = fma(mpole[8], E[ti(l+1,m,n+1,fieldOrder)], total);
  total = fma(mpole[9], E[ti(l,m+1,n+1,fieldOrder)], total);
  return total;
}

void multipoleInteraction(System* system, int i, int j) {
  // Generate Ewald Tensor with the distance between i and j
  int tensorOrder = 5; // Order of the tensor
  // Order of the electric field = tensorOrder-2 since qxx requires d^2/dx^2(1/r)
  int fieldOrder = tensorOrder-2;
  REAL r[3];
  for(int k = 0; k < 3; k++) {
    // Axis lengths
    REAL* a = system->boxDim[k];
    REAL axisLen = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    r[k] = pbc(system->X[i*3+k] - system->X[j*3+k], axisLen);
  }
  REAL src[tensorOrder+1];
  for(int i = 0; i < tensorOrder; i++) {
    src[i] = 0.0;
  }
  ewaldSource(src, r, system->ewaldBeta, tensorOrder);
  int tensorTerms = nChooseK(tensorOrder+3, 3);
  REAL tensor[tensorTerms];
  for(int i = 0; i < tensorTerms; i++) {
    tensor[i] = 0.0;
  }
  generateTensor(tensor, r, src, tensorOrder);

  // Interact multipole I with multipole J
  REAL* mpoleI = system->multipoles[i];
  REAL* mpoleJ = system->multipoles[j];
  REAL E[nChooseK(fieldOrder+3, 3)];
  potentialDueToMultipole(E, fieldOrder, mpoleI, true, tensor, tensorOrder);
  REAL potential = multipoleEnergy(mpoleJ, E);
  REAL xjF = multipoleFieldInteraction(mpoleJ, E, 1, 0, 0, fieldOrder);
  REAL yjF = multipoleFieldInteraction(mpoleJ, E, 0, 1, 0, fieldOrder);
  REAL zjF = multipoleFieldInteraction(mpoleJ, E, 0, 0, 1, fieldOrder);
}
