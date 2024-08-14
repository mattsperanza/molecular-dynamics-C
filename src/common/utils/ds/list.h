#ifndef MOLECULAR_DYNAMICS_C_LIST_H
#define MOLECULAR_DYNAMICS_C_LIST_H

#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>

/**
 * Generic double-linked list using void pointers.
 * <p>
 * C++ creator said these are slow af, so only use if absolutely necessary.
 */
typedef struct _listNode {
    void* data;
    struct _listNode* next;
    struct _listNode* previous;
} ListNode;

typedef void(*CallbackFree)(void *);
typedef int(*CallbackCompare)(void *a, void *b);
typedef bool(*CallbackIterate)(int index, ListNode *node);

typedef struct _list {
    int size;
    int bytesPerData;
    ListNode* head;
    ListNode* tail;
    CallbackFree callbackFree;
    CallbackCompare callbackCompare;
} LinkedList;

// Construct & Destruct
LinkedList* listCreate(int dataSize, CallbackFree free_callback, CallbackCompare compare_callback);
void listDestroy(LinkedList *list);

// Modifiers
void listAdd(LinkedList *list, void *data);
void listRemove(LinkedList *list, int index);

// Search
bool listContains(LinkedList *list, void *data);
ListNode* listFind(LinkedList *list, void *data);
ListNode* listGetNode(LinkedList *list, int index);

// Extraction
void* listGetNodeData(LinkedList *list, int index);
void** listToArray(LinkedList *list, void** ret);

// Misc.
void list_iterate(LinkedList *list, CallbackIterate iterate_callback);
#endif //MOLECULAR_DYNAMICS_C_LIST_H
