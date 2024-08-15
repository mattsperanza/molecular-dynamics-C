#include "../../include/linkedList.h"
#include <assert.h>
#include <stdlib.h>

struct ListNode;
struct LinkedList;

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
