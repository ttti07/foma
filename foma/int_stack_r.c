/*   Foma: a finite-state toolkit and library.                                 */
/*   Copyright © 2008-2015 Mans Hulden                                         */

/*   This file is part of foma.                                                */

/*   Licensed under the Apache License, Version 2.0 (the "License");           */
/*   you may not use this file except in compliance with the License.          */
/*   You may obtain a copy of the License at                                   */

/*      http://www.apache.org/licenses/LICENSE-2.0                             */

/*   Unless required by applicable law or agreed to in writing, software       */
/*   distributed under the License is distributed on an "AS IS" BASIS,         */
/*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/*   See the License for the specific language governing permissions and       */
/*   limitations under the License.                                            */

#include <stdio.h>
#include <stdlib.h>
#include "foma/foma_r.h"

#define MAX_STACK 2097152
#define MAX_PTR_STACK 2097152

struct ptr_stack_handle *ptr_stack_handle_create(int stack_size, int ptr_stack_size) {
    struct ptr_stack_handle *handle = xxmalloc(sizeof(struct ptr_stack_handle));
    if (stack_size <= 0) {
        handle->a = xxmalloc(MAX_STACK * sizeof(int));
        handle->stack_size = MAX_STACK;
    } else {
        handle->a = xxmalloc(stack_size * sizeof(int));
        handle->stack_size = stack_size;
    }
    handle->top = -1;

    if (ptr_stack_size <= 0) {
        handle->ptr_stack = xxmalloc(MAX_PTR_STACK * sizeof(void *));
        handle->ptr_stack_size = MAX_PTR_STACK;
    } else {
        handle->ptr_stack = xxmalloc(ptr_stack_size * sizeof(void *));
        handle->ptr_stack_size = ptr_stack_size;
    }
    handle->ptr_stack_top = -1;

    return handle;
}

void ptr_stack_handle_destroy(struct ptr_stack_handle *handle) {
    xxfree(handle->a);
    xxfree(handle->ptr_stack);
    xxfree(handle);
}

int ptr_stack_isempty(struct ptr_stack_handle *handle) {
    return handle->ptr_stack_top == -1;
}

void ptr_stack_clear(struct ptr_stack_handle *handle) {
    handle->ptr_stack_top = -1;
}

void *ptr_stack_pop(struct ptr_stack_handle *handle) {
    return handle->ptr_stack[handle->ptr_stack_top--];
}

int ptr_stack_isfull(struct ptr_stack_handle *handle) {
    return (handle->ptr_stack_top == (handle->ptr_stack_size - 1));
}

void ptr_stack_push(struct ptr_stack_handle *handle, void *ptr) {
    if (ptr_stack_isfull(handle)) {
        fprintf(stderr, "Pointer stack full!\n");
        exit(1);
    }
    handle->ptr_stack[++(handle->ptr_stack_top)] = ptr;
}


int int_stack_isempty(struct ptr_stack_handle *handle) {
  return handle->top == -1;
}

void int_stack_clear(struct ptr_stack_handle *handle) {
  handle->top = -1;
}

int int_stack_find(struct ptr_stack_handle *handle, int entry) {
  int i;
  if (int_stack_isempty(handle)) {
    return 0;
  }
  for(i = 0; i <= handle->top ; i++) {
    if (entry == handle->a[i]) {
      return 1;
    }
  }
  return 0;
}

int int_stack_size(struct ptr_stack_handle *handle) {
  return (handle->top + 1);
}

void int_stack_push(struct ptr_stack_handle *handle, int c) {
  if (int_stack_isfull(handle)) {
    fprintf(stderr, "Stack full!\n");
    exit(1);
  }
  handle->a[++(handle->top)] = c;
}


int int_stack_pop(struct ptr_stack_handle *handle) {
  return handle->a[handle->top--];
}

int int_stack_isfull(struct ptr_stack_handle *handle) {
  return (handle->top == (handle->stack_size - 1));
}

