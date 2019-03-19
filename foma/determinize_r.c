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
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include "foma/foma_r.h"

#define SUBSET_EPSILON_REMOVE 1
#define SUBSET_DETERMINIZE 2
#define SUBSET_TEST_STAR_FREE 3

#define NHASH_LOAD_LIMIT 2 /* load limit for nhash table size */

struct e_closure_memo {
    int state;
    int mark;
    struct e_closure_memo *target;
    struct e_closure_memo *next;
};

static unsigned int primes[26] = {61,127,251,509,1021,2039,4093,8191,16381,32749,65521,131071,262139,524287,1048573,2097143,4194301,8388593,16777213,33554393,67108859,134217689,268435399,536870909,1073741789,2147483647};

struct nhash_list {
    int setnum;
    unsigned int size;
    unsigned int set_offset;
    struct nhash_list *next;
};

struct T_memo {
    unsigned char finalstart;
    unsigned int size;
    unsigned int set_offset;
};

struct trans_list {
    int inout;
    int target;
};

struct trans_array {
    struct trans_list *transitions;
    unsigned int size;
    unsigned int tail;
};

struct determinize_handle {
    int fsm_linecount, num_states, num_symbols, epsilon_symbol, *single_sigma_array, *double_sigma_array, limit, num_start_states, op;
    _Bool *finals, deterministic, numss;
    struct e_closure_memo *e_closure_memo;
    int T_last_unmarked, T_limit;
    struct trans_list *trans_list;
    struct trans_array *trans_array;
    struct T_memo *T_ptr;
    int nhash_tablesize, nhash_load, current_setnum, *e_table, *marktable, *temp_move, mainloop, maxsigma, *set_table, set_table_size, star_free_mark;
    unsigned int set_table_offset;
    struct nhash_list *table;
    struct build_handle *b_handle;
};

extern int add_fsm_arc(struct fsm_state *fsm, int offset, int state_no, int in, int out, int target, int final_state, int start_state);

static void init(struct determinize_handle *det_handle, struct fsm *net);
INLINE static int e_closure(struct determinize_handle *det_handle, int states);
INLINE static int set_lookup(struct determinize_handle *det_handle, int *lookup_table, int size);
static int initial_e_closure(struct determinize_handle *det_handle, struct fsm *network);
static void memoize_e_closure(struct determinize_handle *det_handle, struct fsm_state *fsm);
static int next_unmarked(struct determinize_handle *det_handle);
static void single_symbol_to_symbol_pair(struct determinize_handle *det_handle, int symbol, int *symbol_in, int *symbol_out);
static int symbol_pair_to_single_symbol(struct determinize_handle *det_handle, int in, int out);
static void sigma_to_pairs(struct determinize_handle *det_handle, struct fsm *net);
static int nhash_find_insert(struct determinize_handle *det_handle, int *set, int setsize);
INLINE static int hashf(struct determinize_handle *det_handle, int *set, int setsize);
static int nhash_insert(struct determinize_handle *det_handle, int hashval, int *set, int setsize);
static void nhash_rebuild_table (struct determinize_handle *det_handle);
static void nhash_init (struct determinize_handle *det_handle, int initial_size);
static void nhash_free(struct nhash_list *nptr, int size);
static void e_closure_free(struct determinize_handle *det_handle);
static void init_trans_array(struct determinize_handle *det_handle, struct fsm *net);
static struct fsm *fsm_subset(struct determinize_handle *det_handle, struct fsm *net, int operation);

struct fsm *fsm_epsilon_remove(struct build_handle *b_handle, struct fsm *net) {
    struct determinize_handle det_handle;
    det_handle.b_handle = b_handle;

    return(fsm_subset(&det_handle, net, SUBSET_EPSILON_REMOVE));
}

struct fsm *fsm_determinize(struct build_handle *b_handle, struct fsm *net) {
    struct determinize_handle det_handle;
    det_handle.b_handle = b_handle;

    return(fsm_subset(&det_handle, net, SUBSET_DETERMINIZE));
}

int fsm_isstarfree(struct fsm *net) {
    #define DFS_WHITE 0
    #define DFS_GRAY 1
    #define DFS_BLACK 2

    struct fsm *sfnet;
    struct state_array *state_array;
    struct fsm_state *curr_ptr;
    int v, vp, is_star_free;
    short int in;
    char *dfs_map;

    struct determinize_handle det_handle;
    sfnet = fsm_subset(&det_handle, net, SUBSET_TEST_STAR_FREE);
    is_star_free = 1;

    state_array = map_firstlines(net);
    ptr_stack_clear(det_handle.b_handle->stack);
    ptr_stack_push(det_handle.b_handle->stack, state_array->transitions);

    dfs_map = xxcalloc(sfnet->statecount, sizeof(char));
    while(!ptr_stack_isempty(det_handle.b_handle->stack)) {

        curr_ptr = ptr_stack_pop(det_handle.b_handle->stack);
    nopop:
        v = curr_ptr->state_no; /* source state number */
        vp = curr_ptr->target;  /* target state number */

        if (v == -1 || vp == -1) {
            *(dfs_map+v) = DFS_BLACK;
            continue;
        }
        *(dfs_map+v) = DFS_GRAY;

        in = curr_ptr->in;
        if (*(dfs_map+vp) == DFS_GRAY && in == det_handle.maxsigma) {
            /* Not star-free */
            is_star_free = 0;
            break;
        }
        if (v == (curr_ptr+1)->state_no) {
            ptr_stack_push(det_handle.b_handle->stack, curr_ptr+1);
        }
        if (*(dfs_map+vp) == DFS_WHITE) { 
            curr_ptr = (state_array+vp)->transitions;
            goto nopop;
        }
    }
    ptr_stack_clear(det_handle.b_handle->stack);
    xxfree(dfs_map);
    xxfree(state_array);
    //stack_add(sfnet);
    return(is_star_free);
}

static struct fsm *fsm_subset(struct determinize_handle *det_handle, struct fsm *net, int operation) {

    int T, U;
    
    if (net->is_deterministic == YES && operation != SUBSET_TEST_STAR_FREE) {
        return(net);
    }
    /* Export this var */
    det_handle->op = operation;
    fsm_count(net);
    det_handle->num_states = net->statecount;
    det_handle->deterministic = 1;
    init(det_handle, net);
    nhash_init(det_handle, (det_handle->num_states < 12) ? 6 : det_handle->num_states/2);
    
    T = initial_e_closure(det_handle, net);

    int_stack_clear(det_handle->b_handle->stack);
    
    if (det_handle->deterministic == 1 && det_handle->epsilon_symbol == -1 && det_handle->num_start_states == 1 && det_handle->numss == 0) {
        net->is_deterministic = YES;
        net->is_epsilon_free = YES;
        nhash_free(det_handle->table, det_handle->nhash_tablesize);
        xxfree(det_handle->T_ptr);
        xxfree(det_handle->e_table);
        xxfree(det_handle->trans_list);
        xxfree(det_handle->trans_array);
        xxfree(det_handle->double_sigma_array);
        xxfree(det_handle->single_sigma_array);
        xxfree(det_handle->finals);
        xxfree(det_handle->temp_move);
        xxfree(det_handle->set_table);
        return(net);
    }

    if (operation == SUBSET_EPSILON_REMOVE && det_handle->epsilon_symbol == -1) {
        net->is_epsilon_free = YES;
        nhash_free(det_handle->table, det_handle->nhash_tablesize);
        xxfree(det_handle->T_ptr);
        xxfree(det_handle->e_table);
        xxfree(det_handle->trans_list);
        xxfree(det_handle->trans_array);
        xxfree(det_handle->double_sigma_array);
        xxfree(det_handle->single_sigma_array);
        xxfree(det_handle->finals);
        xxfree(det_handle->temp_move);
        xxfree(det_handle->set_table);
        return(net);
    }

    if (operation == SUBSET_TEST_STAR_FREE) {
        fsm_state_init(det_handle->b_handle->da_handle, sigma_max(&net->sigma)+1);
        det_handle->star_free_mark = 0;
    } else {
        fsm_state_init(det_handle->b_handle->da_handle, sigma_max(&net->sigma));
        xxfree(net->states);
    }

    /* init */

    do {
        int i, j, tail, setsize, *theset, stateno, has_trans, minsym, next_minsym, trgt, symbol_in, symbol_out;
        struct trans_list *transitions;
        struct trans_array *tptr;

        fsm_state_set_current_state(det_handle->b_handle->da_handle, T, (det_handle->T_ptr+T)->finalstart, T == 0 ? 1 : 0);
        
        /* Prepare set */
        setsize = (det_handle->T_ptr+T)->size;
        theset = det_handle->set_table+(det_handle->T_ptr+T)->set_offset;
        minsym = INT_MAX;
        has_trans = 0;
        for (i = 0; i < setsize; i++) {
            stateno = *(theset+i);
            tptr = det_handle->trans_array+stateno;
            tptr->tail = 0;
            if (tptr->size == 0)
                continue;
            if ((tptr->transitions)->inout < minsym) {
                minsym = (tptr->transitions)->inout;
                has_trans = 1;
            }
        }
        if (!has_trans) {
            /* close state */
            fsm_state_end_state(det_handle->b_handle->da_handle);
            continue;
        }
        
        /* While set not empty */

        for (next_minsym = INT_MAX; minsym != INT_MAX ; minsym = next_minsym, next_minsym = INT_MAX) {
            theset = det_handle->set_table+(det_handle->T_ptr+T)->set_offset;
            
            for (i = 0, j = 0 ; i < setsize; i++) {
                
                stateno = *(theset+i);
                tptr = det_handle->trans_array+stateno;
                tail = tptr->tail;
                transitions = (tptr->transitions)+tail;
                
                while (tail < tptr->size &&  transitions->inout == minsym) {
                    trgt = transitions->target;
                    if (*(det_handle->e_table+(trgt)) != det_handle->mainloop) {
                        *(det_handle->e_table+(trgt)) = det_handle->mainloop;
                        *(det_handle->temp_move+j) = trgt;
                        j++;
                        
                        if (operation == SUBSET_EPSILON_REMOVE) {
                            det_handle->mainloop++;
                            if ((U = e_closure(det_handle, j)) != -1) {
                                single_symbol_to_symbol_pair(det_handle, minsym, &symbol_in, &symbol_out);
                                fsm_state_add_arc(det_handle->b_handle->da_handle, T, symbol_in, symbol_out, U, (det_handle->T_ptr+T)->finalstart, T == 0 ? 1 : 0);
                                j = 0;
                            }
                        }
                    }
                    tail++;
                    transitions++;
                }
                
                tptr->tail = tail;
                
                if (tail == tptr->size)
                    continue;
                /* Check next minsym */
                if (transitions->inout < next_minsym) {
                    next_minsym = transitions->inout;
                }
            }
            if (operation == SUBSET_DETERMINIZE) {
                det_handle->mainloop++;
                if ((U = e_closure(det_handle, j)) != -1) {
                    single_symbol_to_symbol_pair(det_handle, minsym, &symbol_in, &symbol_out);
                    fsm_state_add_arc(det_handle->b_handle->da_handle, T, symbol_in, symbol_out, U, (det_handle->T_ptr+T)->finalstart, T == 0 ? 1 : 0);
                }
            }
            if (operation == SUBSET_TEST_STAR_FREE) {
                det_handle->mainloop++;
                if ((U = e_closure(det_handle, j)) != -1) {
                    single_symbol_to_symbol_pair(det_handle, minsym, &symbol_in, &symbol_out);                   
                    fsm_state_add_arc(det_handle->b_handle->da_handle, T, symbol_in, symbol_out, U, (det_handle->T_ptr+T)->finalstart, T == 0 ? 1 : 0);
                    if (det_handle->star_free_mark == 1) {
                        //fsm_state_add_arc(det_handle->b_handle->da_handle, T, det_handle->maxsigma, det_handle->maxsigma, U, (det_handle->T_ptr+T)->finalstart, T == 0 ? 1 : 0);
                        det_handle->star_free_mark = 0;
                    }
                }
            }
        }
        /* end state */
        fsm_state_end_state(det_handle->b_handle->da_handle);
    } while ((T = next_unmarked(det_handle)) != -1);
    
    /* wrapup() */
    nhash_free(det_handle->table, det_handle->nhash_tablesize);
    xxfree(det_handle->set_table);
    xxfree(det_handle->T_ptr);
    xxfree(det_handle->temp_move);
    xxfree(det_handle->e_table);
    xxfree(det_handle->trans_list);
    xxfree(det_handle->trans_array);
    
    if (det_handle->epsilon_symbol != -1)
        e_closure_free(det_handle);
    xxfree(det_handle->double_sigma_array);
    xxfree(det_handle->single_sigma_array);
    xxfree(det_handle->finals);
    fsm_state_close(det_handle->b_handle->da_handle, net);
    return(net);
}

static void init(struct determinize_handle *det_handle, struct fsm *net) {
    /* A temporary table for handling epsilon closure */
    /* to avoid doubles */

    det_handle->e_table = xxcalloc(net->statecount,sizeof(int));

    /* Counter for our access tables */

    det_handle->mainloop = 1;

    /* Temporary table for storing sets and */
    /* passing to hash function */
    
    /* Table for listing current results of move & e-closure */
    det_handle->temp_move = xxmalloc((net->statecount + 1) *sizeof(int));
    
    /* We malloc this much memory to begin with for the new fsm */
    /* Then grow it by the double as needed */

    det_handle->limit = next_power_of_two(net->linecount);
    det_handle->fsm_linecount = 0;
    sigma_to_pairs(det_handle, net);

    /* Optimistically malloc det_handle->T_ptr array */
    /* We allocate memory for a number of pointers to a set of states */
    /* To handle fast lookup in array */
    /* Optimistically, we choose the initial size to be the number of */
    /* states in the non-deterministic fsm */
    
    det_handle->T_last_unmarked = 0;
    det_handle->T_limit = next_power_of_two(det_handle->num_states);

    det_handle->T_ptr = xxcalloc(det_handle->T_limit,sizeof(struct T_memo));

    /* Stores all sets consecutively in one table */
    /* det_handle->T_ptr->set_offset and size                 */
    /* are used to retrieve the set               */

    det_handle->set_table_size = next_power_of_two(det_handle->num_states);
    det_handle->set_table = xxmalloc(det_handle->set_table_size*sizeof(int));
    det_handle->set_table_offset = 0;

    init_trans_array(det_handle, net);
}

static int trans_sort_cmp(const void *a, const void *b) {
  return (((const struct trans_list *)a)->inout - ((const struct trans_list *)b)->inout);
}

static void init_trans_array(struct determinize_handle *det_handle, struct fsm *net) {
    struct trans_list *arrptr;
    struct fsm_state *fsm;
    int i, j, laststate, lastsym, inout, size, state;

    arrptr = det_handle->trans_list = xxmalloc(net->linecount * sizeof(struct trans_list));
    det_handle->trans_array = xxcalloc(net->statecount, sizeof(struct trans_array));
    
    laststate = -1;
    fsm = net->states;

    for (i=0, size = 0; (fsm+i)->state_no != -1; i++) {
        state = (fsm+i)->state_no;
        if (state != laststate) {
            if (laststate != -1) {
                (det_handle->trans_array+laststate)->size = size;
            }
            (det_handle->trans_array+state)->transitions = arrptr;
            size = 0;
        }
        laststate = state;

        if ((fsm+i)->target == -1)
            continue;
        inout = symbol_pair_to_single_symbol(det_handle, (fsm+i)->in, (fsm+i)->out);
        if (inout == det_handle->epsilon_symbol)
            continue;
        
        arrptr->inout = inout;
        arrptr->target = (fsm+i)->target;
        arrptr++;
        size++;
    }

    if (laststate != -1) {
        (det_handle->trans_array+laststate)->size = size;
    }

    for (i=0; i < net->statecount; i++) {
        arrptr = (det_handle->trans_array+i)->transitions;
        size = (det_handle->trans_array+i)->size;
        if (size > 1) {
            qsort(arrptr, size, sizeof(struct trans_list), trans_sort_cmp);
            lastsym = -1;
            /* Figure out if we're already deterministic */
            for (j=0; j < size; j++) {
                if ((arrptr+j)->inout == lastsym)
                    det_handle->deterministic = 0;
                lastsym = (arrptr+j)->inout;
            }
        }
    }
}

INLINE static int e_closure(struct determinize_handle *det_handle, int states) {

    int i, set_size;
    struct e_closure_memo *ptr;

    /* e_closure extends the list of states which are reachable */
    /* and appends these to det_handle->e_table                             */

    if (det_handle->epsilon_symbol == -1)
        return(set_lookup(det_handle, det_handle->temp_move, states));

    if (states == 0)
        return -1;

    det_handle->mainloop--;
    
    set_size = states;

    for (i = 0; i < states; i++) {

        /* State number we want to do e-closure on */
        ptr = det_handle->e_closure_memo + *(det_handle->temp_move+i);
        if (ptr->target == NULL)
            continue;
        ptr_stack_push(det_handle->b_handle->stack, ptr);

        while (!(ptr_stack_isempty(det_handle->b_handle->stack))) {
            ptr = ptr_stack_pop(det_handle->b_handle->stack);
            /* Don't follow if already seen */
            if (*(det_handle->marktable+ptr->state) == det_handle->mainloop)
                continue;
            
            ptr->mark = det_handle->mainloop;
            *(det_handle->marktable+ptr->state) = det_handle->mainloop;
            /* Add to tail of list */
            if (*(det_handle->e_table+(ptr->state)) != det_handle->mainloop) {
                *(det_handle->temp_move+set_size) = ptr->state;
                *(det_handle->e_table+(ptr->state)) = det_handle->mainloop;
                set_size++;
            }
            
            if (ptr->target == NULL)
                continue;
            /* Traverse chain */

            for (; ptr != NULL ; ptr = ptr->next) {
                if (ptr->target->mark != det_handle->mainloop) {
                    /* Push */
                    ptr->target->mark = det_handle->mainloop;
                    ptr_stack_push(det_handle->b_handle->stack, ptr->target);
                }
            }
        }
    }

    det_handle->mainloop++;
    return(set_lookup(det_handle, det_handle->temp_move, set_size));
}

INLINE static int set_lookup (struct determinize_handle *det_handle, int *lookup_table, int size) {

  /* Look up a set and its corresponding state number */
  /* if it doesn't exist from before, assign a state number */
  
    return(nhash_find_insert(det_handle, lookup_table, size));
  
}

static void add_T_ptr(struct determinize_handle *det_handle, int setnum, int setsize, unsigned int theset, int fs) {

  int i;
  if (setnum >= det_handle->T_limit) {
    det_handle->T_limit *= 2;
    det_handle->T_ptr = xxrealloc(det_handle->T_ptr, sizeof(struct T_memo)*det_handle->T_limit);
    for (i=setnum; i < det_handle->T_limit; i++) {
        (det_handle->T_ptr+i)->size = 0;
    }
  }
  
  (det_handle->T_ptr + setnum)->size = setsize;
  (det_handle->T_ptr + setnum)->set_offset = theset;
  (det_handle->T_ptr + setnum)->finalstart = fs;
  int_stack_push(det_handle->b_handle->stack, setnum);

}

static int initial_e_closure(struct determinize_handle *det_handle, struct fsm *net) {

    struct fsm_state *fsm;
    int i,j;

    det_handle->finals = xxcalloc(det_handle->num_states, sizeof(_Bool));

    det_handle->num_start_states = 0;
    fsm = net->states;
    
    /* Create lookups for each state */
    for (i=0,j=0; (fsm+i)->state_no != -1; i++) {
        if ((fsm+i)->final_state) {
            det_handle->finals[(fsm+i)->state_no] = 1;
        }
        /* Add the start states as the initial set */
        if ((det_handle->op == SUBSET_TEST_STAR_FREE) || ((fsm+i)->start_state)) {
            if (*(det_handle->e_table+((fsm+i)->state_no)) != det_handle->mainloop) {
                det_handle->num_start_states++;
                det_handle->numss = (fsm+i)->state_no;
                *(det_handle->e_table+((fsm+i)->state_no)) = det_handle->mainloop;
                *(det_handle->temp_move+j) = (fsm+i)->state_no;
                j++;
            }
        }
    }
    det_handle->mainloop++;
    /* Memoize e-closure(u) */
    if (det_handle->epsilon_symbol != -1) {
        memoize_e_closure(det_handle, fsm);
    }
    return(e_closure(det_handle, j));
}
 
static void memoize_e_closure(struct determinize_handle *det_handle, struct fsm_state *fsm) {
    
    int i, state, laststate, *redcheck;
    struct e_closure_memo *ptr;
    
    det_handle->e_closure_memo = xxcalloc(det_handle->num_states,sizeof(struct e_closure_memo));
    det_handle->marktable = xxcalloc(det_handle->num_states,sizeof(int));
    /* Table for avoiding redundant epsilon arcs in closure */
    redcheck = xxmalloc(det_handle->num_states*sizeof(int));

    for (i=0; i < det_handle->num_states; i++) {
        ptr = det_handle->e_closure_memo+i;
        ptr->state = i;
        ptr->target = NULL;
        *(redcheck+i) = -1;
    }

    laststate = -1;

    for (i=0; ;i++) {
        
        state = (fsm+i)->state_no;
        
        if (state != laststate) {
            if (!int_stack_isempty(det_handle->b_handle->stack)) {                
                det_handle->deterministic = 0;
                ptr = det_handle->e_closure_memo+laststate;
                ptr->target = det_handle->e_closure_memo+int_stack_pop(det_handle->b_handle->stack);
                while (!int_stack_isempty(det_handle->b_handle->stack)) {
                    ptr->next = xxmalloc(sizeof(struct e_closure_memo));
                    ptr->next->state = laststate;
                    ptr->next->target = det_handle->e_closure_memo+int_stack_pop(det_handle->b_handle->stack);
                    ptr->next->next = NULL;
                    ptr = ptr->next;
                }
            }
        }
        if (state == -1) {
            break;
        }
        if ((fsm+i)->target == -1) {
            continue;
        }
        /* Check if we have a redundant epsilon arc */
        if ((fsm+i)->in == EPSILON && (fsm+i)->out == EPSILON) {
            if (*(redcheck+((fsm+i)->target)) != (fsm+i)->state_no) {
                if ((fsm+i)->target != (fsm+i)->state_no) {
                    int_stack_push(det_handle->b_handle->stack, (fsm+i)->target);
                    *(redcheck+((fsm+i)->target)) = (fsm+i)->state_no;
                }
            }
            laststate = state;
        }
    }
    xxfree(redcheck);
}
 
static int next_unmarked(struct determinize_handle *det_handle) {
    if ((int_stack_isempty(det_handle->b_handle->stack)))
        return -1;
    return(int_stack_pop(det_handle->b_handle->stack));

    if ((det_handle->T_limit <= det_handle->T_last_unmarked + 1) || (det_handle->T_ptr+det_handle->T_last_unmarked+1)->size == 0) {
        return -1;
    } else {
        det_handle->T_last_unmarked++;
        return(det_handle->T_last_unmarked);
    }
}

static void single_symbol_to_symbol_pair(struct determinize_handle *det_handle, int symbol, int *symbol_in, int *symbol_out) {

  *symbol_in = *(det_handle->single_sigma_array+(symbol*2));
  *symbol_out = *(det_handle->single_sigma_array+(symbol*2+1));
  
}

static int symbol_pair_to_single_symbol(struct determinize_handle *det_handle, int in, int out) {
  return(*(det_handle->double_sigma_array+det_handle->maxsigma*in+out));
}

static void sigma_to_pairs(struct determinize_handle *det_handle, struct fsm *net) {
  
  int i, j, x, y, z, next_x = 0;
  struct fsm_state *fsm;

  fsm = net->states;

  det_handle->epsilon_symbol = -1; 
  det_handle->maxsigma = sigma_max(&net->sigma);
  det_handle->maxsigma++;

  det_handle->single_sigma_array = xxmalloc(2 * det_handle->maxsigma * det_handle->maxsigma * sizeof(int));
  det_handle->double_sigma_array = xxmalloc(det_handle->maxsigma * det_handle->maxsigma * sizeof(int));
  
  for (i=0; i < det_handle->maxsigma; i++) {
    for (j=0; j< det_handle->maxsigma; j++) {
      *(det_handle->double_sigma_array+det_handle->maxsigma*i+j) = -1;
    }
  }
  
  /* f(x) -> y,z sigma pair */
  /* f(y,z) -> x simple entry */
  /* if exists f(n) <-> EPSILON, EPSILON, save n */
  /* symbol(x) x>=1 */

  /* Forward mapping: */
  /* *(det_handle->double_sigma_array+det_handle->maxsigma*in+out) */

  /* Backmapping: */
  /* *(det_handle->single_sigma_array+(symbol*2) = in(symbol) */
  /* *(det_handle->single_sigma_array+(symbol*2+1) = out(symbol) */

  /* Table for checking whether a state is final */

  x = 0;
  net->arity = 1;
  for (i=0; (fsm+i)->state_no != -1; i++) {
    y = (fsm+i)->in;
    z = (fsm+i)->out;
    if ((y == -1) || (z == -1))
      continue;
    if (y != z || y == UNKNOWN || z == UNKNOWN)
        net->arity = 2;
    if (*(det_handle->double_sigma_array+det_handle->maxsigma*y+z) == -1) {
      *(det_handle->double_sigma_array+det_handle->maxsigma*y+z) = x;
      *(det_handle->single_sigma_array+next_x) = y;
      next_x++;
      *(det_handle->single_sigma_array+next_x) = z;
      next_x++;
      if (y == EPSILON && z == EPSILON) {
	det_handle->epsilon_symbol = x;
      }
      x++;
    }
  }
  det_handle->num_symbols = x;
}


/* Functions for hashing n integers */
/* with permutations hashing to the same value */
/* necessary for subset construction */

static int nhash_find_insert(struct determinize_handle *det_handle, int *set, int setsize) {
    int j, found, *currlist;
    struct nhash_list *tableptr;
    unsigned int hashval;
    
    hashval = hashf(det_handle, set, setsize);
    if ((det_handle->table+hashval)->size == 0) {
        return(nhash_insert(det_handle, hashval, set, setsize));
    } else {
        for (tableptr=(det_handle->table+hashval); tableptr != NULL; tableptr = tableptr->next) {
            if ((tableptr)->size != setsize) {
                continue;
            } else {
                /* Compare the list at hashval position */
                /* to the current set by looking at etable */
                /* entries */
                for (j=0, found = 1, currlist= det_handle->set_table+tableptr->set_offset; j < setsize; j++) {
                    if (*(det_handle->e_table+(*(currlist+j))) != (det_handle->mainloop-1)) {
                        found = 0;
                        break;
                    }
                }
                if (det_handle->op == SUBSET_TEST_STAR_FREE && found == 1) {
                    for (j=0, currlist= det_handle->set_table+tableptr->set_offset; j < setsize; j++) {
                        if (*(set+j) != *(currlist+j)) {
                            /* Set mark */
                            det_handle->star_free_mark = 1;
                        }
                    }
                }
                if (found == 1) {
                    return(tableptr->setnum);
                }
            }
        }
        
        if (det_handle->nhash_load / NHASH_LOAD_LIMIT > det_handle->nhash_tablesize) {
            nhash_rebuild_table(det_handle);
            hashval = hashf(det_handle, set, setsize);
        }
        return(nhash_insert(det_handle, hashval, set, setsize));
    }
}

INLINE static int hashf(struct determinize_handle *det_handle, int *set, int setsize) {
  int i;
  unsigned int hashval, sum = 0;
  hashval = 6703271;
  for (i = 0; i < setsize; i++) {
      hashval = (unsigned int) (*(set+i) + 1103 * setsize) * hashval; 
      sum += *(set+i) + i;
  }
  hashval = hashval + sum * 31;
  hashval = (hashval % det_handle->nhash_tablesize);
  return hashval;
}

static unsigned int move_set(struct determinize_handle *det_handle, int *set, int setsize) {
    unsigned int old_offset;
    if (det_handle->set_table_offset + setsize >= det_handle->set_table_size) {
        while (det_handle->set_table_offset + setsize >= det_handle->set_table_size) {
            det_handle->set_table_size *= 2;
        }
        det_handle->set_table = xxrealloc(det_handle->set_table, det_handle->set_table_size * sizeof(int));
    }
    memcpy(det_handle->set_table+det_handle->set_table_offset, set, setsize * sizeof(int));
    old_offset = det_handle->set_table_offset;
    det_handle->set_table_offset += setsize;
    return(old_offset);
}

static int nhash_insert(struct determinize_handle *det_handle, int hashval, int *set, int setsize) { 
  struct nhash_list *tableptr;  
  int i, fs = 0;

  det_handle->current_setnum++;
  tableptr = det_handle->table+hashval;

  det_handle->nhash_load++;
  for (i = 0; i < setsize; i++) {
      if (det_handle->finals[*(set+i)])
          fs = 1;
  }
  if (tableptr->size == 0) {
    
      tableptr->set_offset = move_set(det_handle, set, setsize);
      tableptr->size = setsize;
      tableptr->setnum = det_handle->current_setnum;
      
      add_T_ptr(det_handle, det_handle->current_setnum, setsize, tableptr->set_offset, fs);
      return(det_handle->current_setnum);
  }
  
  tableptr = xxmalloc(sizeof(struct nhash_list));
  tableptr->next = (det_handle->table+hashval)->next;
  (det_handle->table+hashval)->next = tableptr;
  tableptr->setnum = det_handle->current_setnum;
  tableptr->size = setsize;
  tableptr->set_offset = move_set(det_handle, set, setsize);
  
  add_T_ptr(det_handle, det_handle->current_setnum, setsize, tableptr->set_offset, fs);
  return(det_handle->current_setnum);
}

static void nhash_rebuild_table (struct determinize_handle *det_handle) {
    int i, oldsize;
    struct nhash_list *oldtable, *tableptr, *ntableptr, *newptr;
    unsigned int hashval;
    
    oldtable = det_handle->table;
    oldsize = det_handle->nhash_tablesize;

    det_handle->nhash_load = 0;
    for (i=0; primes[i] < det_handle->nhash_tablesize; i++) { }
    det_handle->nhash_tablesize = primes[(i+1)];
    
    det_handle->table = xxcalloc(det_handle->nhash_tablesize,sizeof(struct nhash_list));
    for (i=0; i < oldsize;i++) {
        if ((oldtable+i)->size == 0) {
            continue;
        }
        tableptr = oldtable+i;
        for ( ; tableptr != NULL; (tableptr = tableptr->next)) {
            /* rehash */
            hashval = hashf(det_handle, det_handle->set_table+tableptr->set_offset,tableptr->size);
            ntableptr = det_handle->table+hashval;
            if (ntableptr->size == 0) {
                det_handle->nhash_load++;
                ntableptr->size = tableptr->size;
                ntableptr->set_offset = tableptr->set_offset;
                ntableptr->setnum = tableptr->setnum;
                ntableptr->next = NULL;
            } else {
                newptr = xxmalloc(sizeof(struct nhash_list));
                newptr->next = ntableptr->next;
                ntableptr->next = newptr;
                newptr->setnum = tableptr->setnum;
                newptr->size = tableptr->size;
                newptr->set_offset = tableptr->set_offset;
            }
        }
    }
    nhash_free(oldtable, oldsize);
}

static void nhash_init (struct determinize_handle *det_handle, int initial_size) {

  int i;

  for (i=0; primes[i] < initial_size; i++) { }
  det_handle->nhash_load = 0;
  det_handle->nhash_tablesize = primes[i];
  det_handle->table = xxcalloc(det_handle->nhash_tablesize , sizeof(struct nhash_list));
  det_handle->current_setnum = -1;
}
static void e_closure_free(struct determinize_handle *det_handle) {
    int i;
    struct e_closure_memo *eptr, *eprev;
    xxfree(det_handle->marktable);
    for (i=0;i < det_handle->num_states; i++) {
        eptr = (det_handle->e_closure_memo+i)->next;
        for (eprev = NULL; eptr != NULL; ) {
            eprev = eptr;
            eptr = eptr->next;
            xxfree(eprev);
        }
        
    }
    xxfree(det_handle->e_closure_memo);
}

static void nhash_free(struct nhash_list *nptr, int size) {
    struct nhash_list *nptr2, *nnext;
    int i;
    for (i=0; i < size; i++) {
        for (nptr2 = (nptr+i)->next; nptr2 != NULL; nptr2 = nnext) {
            nnext = nptr2->next;
            xxfree(nptr2);
        }
    }
    xxfree(nptr);
}

