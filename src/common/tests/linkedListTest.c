/**
 * Uses asserts to test code.
 */
#include <stdio.h>
#include "../utils/ds/linkedList.c"

void printIntegerList(LinkedList const* list);
int compareInt(void* intOne, void* intTwo);

void linkedListTest() {
 // LinkedList test with integers
 LinkedList* intList = listCreate(sizeof(int), compareInt, NULL);
 for(int i = 0; i < 10; i++) {
  int* newInt = malloc(sizeof(int));
  *newInt = i;
  listAppend(intList, newInt);
 }
 // TEST: Remove/Add (added above)
 int firstElement = *(int*)intList->head->data;
 assert(firstElement == 0);
 int lastElement = *(int*)intList->tail->data;
 assert(lastElement == 9);

 //printf("Initial list: ");
 //printIntegerList(intList);

 listRemove(intList, 0);
 //printf("Remove first: ");
 //printIntegerList(intList);

 listRemove(intList, 1);
 //printf("Remove index1:");
 //printIntegerList(intList);

 listRemove(intList, intList->size-1);
 //printf("Remove last:  ");
 //printIntegerList(intList);
 assert(intList->size == 7);

 // TEST: listToArray()
 int listSize = intList->size;
 int listAsArrayExpected[7] = {1, 3, 4, 5, 6, 7, 8};
 int** listAsArray = malloc(listSize*sizeof(int*));
 listToArray(intList, (void**) listAsArray, listSize);
 for(int i = 0; i < listSize; i++) {
  int val1 = listAsArrayExpected[i];
  int val2 = *listAsArray[i];
  assert(val1 == val2);
 }

 // TEST: listFind() & listGetNodeData() & listGetNode
 ListNode* node = listFind(intList, listAsArray[3]);
 assert(*(int*)node->data == *(int*)listAsArray[3]);
 node = intList->tail;
 int ret = *(int*)listGetNodeData(intList, listSize-1);
 assert(*(int*)node->data == ret);
 node = listGetNode(intList, listSize-1);
 assert(*(int*)node->data == ret);

 // TEST: listContains()
 assert(listContains(intList, &ret));

 // Clean
 free(listAsArray);
 listDestroy(intList);

 printf("All tests in file src/common/tests/linkedListTest.c passed!");
}

void printIntegerList(LinkedList const* list){
 assert(list != NULL);
 ListNode * node = list->head;
 for(int i = 0; i < list->size; i++) {
  int data = *(int*)node->data;
  printf("%2d->", data);
  node = node->next;
 }
 printf("\n");
}

int compareInt(void* intOne, void* intTwo) {
 int a = *(int*)intOne;
 int b = *(int*)intTwo;
 if(a < b) {
  return -1;
 }
 if(a == b) {
  return 0;
 }
 return 1;
}
