#include "../../include/linkedList.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// Question
bool hasNext(ListNode *node) {
    return !(node == NULL || node->next == NULL);
}

bool hasPrevious(ListNode *node) {
    return !(node == NULL || node->previous == NULL);
}

/**
 * Compares underlying data of the list.
 * @param list list to search
 * @param data data pointer that gets dereferenced
 * @return if the data pointers' data is pointed to in the list, null if compare callback is not provided
 */
bool listContains(LinkedList *list, void *data) {
    assert(list != NULL);
    assert(data != NULL);
    if (list->callbackCompare == NULL) {
        return NULL;
    }
    ListNode *current = list->head;
    while (hasNext(current) || current == list->tail) {
        if (list->callbackCompare(data, current->data) == 0) {
            return true;
        }
        current = (ListNode *) current->next;
    }
    return false;
}

// Modifiers
/**
 * Add node to the end of the list.
 * @param list list to be appended to
 * @param data data to store in new node at the end of the list
 */
void listAppend(LinkedList *list, void *data) {
    assert(list != NULL);
    // New ptr and add data ptr
    ListNode *nodePtr = malloc(sizeof(ListNode));
    nodePtr->data = data;
    // List adjustments
    if (list->size == 0) {
        list->head = nodePtr;
        nodePtr->previous = NULL;
    } else {
        ListNode* currentTail = list->tail;
        currentTail->next = nodePtr;
        nodePtr->previous = currentTail;
    }
    nodePtr->next = NULL;
    list->tail = nodePtr;
    list->size++;
}

/**
 * @param list list to modify
 * @param index element to remove
 * @return whether the element was successfully removed
 */
bool listRemove(LinkedList *list, int index) {
    assert(list != NULL);
    if (index >= list->size) {
        return false;
    }
    // Find node to remove
    ListNode *current = list->head;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    // Unlink adjacent nodes
    ListNode *prev = current->previous;
    ListNode *next = current->next;
    if(prev != NULL && prev->next != NULL) {
        prev->next = next;
    }
    if(next != NULL && next->previous != NULL) {
        next->previous = prev;
    }
    // Fix head/tail -> sets to null if only one node
    if(current == list->head) {
        list->head = next;
    }
    if(current == list->tail) {
        list->tail = prev;
    }
    // Free data then free node struct
    if (list->callbackFree != NULL) {
        list->callbackFree(current->data);
    } else {
        free(current->data);
    }
    free(current);
    list->size--;
    return true;
}

// Search
ListNode* listFind(LinkedList* list, void* data) {
    assert(list != NULL);
    assert(list->callbackCompare != NULL);
    ListNode *current = list->head;
    while (hasNext(current) || current == list->tail) {
        if (list->callbackCompare(data, current->data) == 0) {
            return current;
        }
        current = current->next;
    }
    return NULL;
}

/**
 * @param list list to be searched
 * @param index node index to return
 * @return node at this list index
 */
ListNode *listGetNode(LinkedList *list, int index) {
    assert(list != NULL);
    if (index >= list->size) {
        return NULL;
    }
    ListNode *current = list->head;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    return current;
}

// Extraction
/**
 * @param list list to fetch from
 * @param index node index to retrieve the data
 * @return data from node at that index
 */
void* listGetNodeData(LinkedList *list, int index) {
    assert(list != NULL);
    if (index >= list->size) {
        return NULL;
    }
    ListNode *current = list->head;
    for (int i = 0; i < index; i++) {
        current = current->next;
    }
    return current->data;
}

/**
 * @param list list to extract data from
 * @param ret pre-allocated array with length = size of list
 * @param size size of the list
 */
void listToArray(LinkedList *list, void **ret, int size) {
    assert(list != NULL);
    ListNode const *current = list->head;
    for (int i = 0; i < size; i++) {
        ret[i] = current->data;
        current = current->next;
    }
}

// Construct & Destroy
LinkedList *listCreate(int dataSize, CallbackCompare cbCompare, CallbackFree cbFree) {
    LinkedList *list = malloc(sizeof(LinkedList));
    list->size = 0;
    list->bytesPerData = dataSize;
    list->head = NULL;
    list->tail = NULL;
    list->callbackCompare = cbCompare;
    list->callbackFree = cbFree;
    return list;
};

/**
 * Deletes list, all its nodes, and all the node data.
 * @param list list to destroy
 */
void listDestroy(LinkedList *list) {
    assert(list != NULL);
    ListNode *current = list->head;
    for (int i = 0; i < list->size-1; i++) {
        ListNode *temp = current;
        // shift pointer and use temp to free
        current = current->next;
        // Free data with passed in
        if (list->callbackFree != NULL) {
            list->callbackFree(temp->data);
        } else {
            free(temp->data);
        }
        free(temp);
    }
    free(list);
};

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

//////////////////////////////////////////// TESTS

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

void linkedListTest(bool verbose) {
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
 if(verbose) {
     printf("Initial item list: ");
     printIntegerList(intList);
 }

 listRemove(intList, 0);
 if(verbose) {
     printf("Remove index 0: ");
     printIntegerList(intList);
 }
 listRemove(intList, 1);
 if(verbose) {
     printf("Remove index 1: ");
     printIntegerList(intList);
 }
 listRemove(intList, intList->size-1);
 if(verbose) {
     printf("Remove last id: ");
     printIntegerList(intList);
 }
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

 printf("All tests of linkedList.c passed!\n");
}

