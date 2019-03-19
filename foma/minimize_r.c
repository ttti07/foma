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

#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include "foma/foma_r.h"

struct statesym {
    int target;
    unsigned short int symbol;
    struct state_list *states;
    struct statesym *next;
};

struct state_list {
    int state;
    struct state_list *next;
};

struct p {
    struct e *first_e;
    struct e *last_e;
    struct p *current_split;
    struct p *next;
    struct agenda *agenda;
    int count;
    int t_count;
    int inv_count;
    int inv_t_count;
};

struct e {
  struct p *group;
  struct e *left;
  struct e *right;
  int inv_count;
};

struct agenda {
  struct p *p;
  struct agenda *next;
  _Bool index;
};

struct trans_list {
    int inout;
    int source;
};

struct trans_array {
    struct trans_list *transitions;
    unsigned int size;
    unsigned int tail;
};

struct minimize_handle {
    int *single_sigma_array, *double_sigma_array, *memo_table, *temp_move, *temp_group, maxsigma, epsilon_symbol, num_states, num_symbols, num_finals, mainloop, total_states;
    _Bool *finals;
    struct trans_list *trans_list;
    struct trans_array *trans_array;
    struct p *P, *Phead, *Pnext, *current_w;
    struct e *E;
    struct agenda *Agenda_head, *Agenda_top, *Agenda_next, *Agenda;
};

static struct fsm *fsm_minimize_brz(struct build_handle *b_handle, struct fsm *net);
static struct fsm *fsm_minimize_hop(struct minimize_handle *min_handle, struct fsm *net);
static struct fsm *rebuild_machine(struct minimize_handle *min_handle, struct fsm *net);

static INLINE int refine_states(struct minimize_handle *min_handle, int sym);
static void init_PE(struct minimize_handle *min_handle);
static void agenda_add(struct minimize_handle *min_handle, struct p *pptr, int start);
static void sigma_to_pairs(struct minimize_handle *min_handle, struct fsm *net);
/* static void single_symbol_to_symbol_pair(int symbol, int *symbol_in, int *symbol_out); */
static INLINE int symbol_pair_to_single_symbol(struct minimize_handle *min_handle, int in, int out);
static void generate_inverse(struct minimize_handle *min_handle, struct fsm *net);

struct fsm *fsm_minimize(struct build_handle *b_handle, struct fsm *net) {
    extern int g_minimal;
    extern int g_minimize_hopcroft;

    if (net == NULL) { return NULL; }
    /* The network needs to be deterministic and trim before we minimize */
    if (net->is_deterministic != YES)
        net = fsm_determinize(b_handle, net);
    if (net->is_pruned != YES)
        net = fsm_coaccessible(b_handle, net);
    if (net->is_minimized != YES && g_minimal == 1) {
        if (g_minimize_hopcroft != 0) {
            struct minimize_handle min_handle;
            net = fsm_minimize_hop(&min_handle, net);
        }
        else 
            net = fsm_minimize_brz(b_handle, net);        
        fsm_update_flags(net,YES,YES,YES,YES,UNK,UNK);
    }
    return(net);
}

static struct fsm *fsm_minimize_brz(struct build_handle *b_handle, struct fsm *net) {
    return(fsm_determinize(b_handle, fsm_reverse(b_handle, fsm_determinize(b_handle, fsm_reverse(b_handle, net)))));
}

static struct fsm *fsm_minimize_hop(struct minimize_handle *min_handle, struct fsm *net) {

    struct e *temp_E;
    struct trans_array *tptr;
    struct trans_list *transitions;
    int i,j,minsym,next_minsym,current_i, stateno, thissize, source;  
    unsigned int tail;

    fsm_count(net);
    if (net->finalcount == 0)  {
	fsm_destroy(net);
	return(fsm_empty_set());
    }

    min_handle->num_states = net->statecount;
    
    min_handle->P = NULL;

    /* 
       1. generate the inverse lookup table
       2. generate P and E (partitions, states linked list)
       3. Init Agenda = {Q, Q-F}
       4. Split until Agenda is empty
    */
    
    sigma_to_pairs(min_handle, net);
    
    init_PE(min_handle);

    if (min_handle->total_states == min_handle->num_states) {
        goto bail;
    }

    generate_inverse(min_handle, net);


    min_handle->Agenda_head->index = 0;
    if (min_handle->Agenda_head->next != NULL)
        min_handle->Agenda_head->next->index = 0;

    for (min_handle->Agenda = min_handle->Agenda_head; min_handle->Agenda != NULL; ) {
        /* Remove current_w from agenda */
        min_handle->current_w = min_handle->Agenda->p;
        current_i = min_handle->Agenda->index;
        min_handle->Agenda->p->agenda = NULL;
        min_handle->Agenda = min_handle->Agenda->next;

        /* Store current group state number in tmp_group */
        /* And figure out minsym */
        /* If index is 0 we start splitting from the first symbol */
        /* Otherwise we split from where we left off last time */

        thissize = 0;
        minsym = INT_MAX;
        for (temp_E = min_handle->current_w->first_e; temp_E != NULL; temp_E = temp_E->right) {
            stateno = temp_E - min_handle->E;
            *(min_handle->temp_group+thissize) = stateno;
            thissize++;
            tptr = min_handle->trans_array+stateno;
            /* Clear tails if symloop should start from 0 */
            if (current_i == 0)
                tptr->tail = 0;
            
            tail = tptr->tail;
            transitions = (tptr->transitions)+tail;
            if (tail < tptr->size && transitions->inout < minsym) {
                minsym = transitions->inout;
            }
        }

        for (next_minsym = INT_MAX; minsym != INT_MAX ; minsym = next_minsym, next_minsym = INT_MAX) {

            /* Add states to min_handle->temp_move */
            for (i = 0, j = 0; i < thissize; i++) {
                tptr = min_handle->trans_array+*(min_handle->temp_group+i);
                tail = tptr->tail;
                transitions = (tptr->transitions)+tail;
                while (tail < tptr->size && transitions->inout == minsym) {
                    source = transitions->source;
                    if (*(min_handle->memo_table+(source)) != min_handle->mainloop) {
                        *(min_handle->memo_table+(source)) = min_handle->mainloop;
                        *(min_handle->temp_move+j) = source;
                        j++;
                    }
                    tail++;
                    transitions++;
                }
                tptr->tail = tail;
                if (tail < tptr->size && transitions->inout < next_minsym) {
                    next_minsym = transitions->inout;
                }
            }
            if (j == 0) {
                continue;
            }
            min_handle->mainloop++;
            if (refine_states(min_handle, j) == 1) {
                break; /* break loop if we split current_w */
            }
        }
        if (min_handle->total_states == min_handle->num_states) {
            break;
        }
    }

    net = rebuild_machine(min_handle, net);

    xxfree(min_handle->trans_array);
    xxfree(min_handle->trans_list);

 bail:
    
    xxfree(min_handle->Agenda_top);
    
    xxfree(min_handle->memo_table);
    xxfree(min_handle->temp_move);
    xxfree(min_handle->temp_group);


    xxfree(min_handle->finals);
    xxfree(min_handle->E);
    xxfree(min_handle->Phead);
    xxfree(min_handle->single_sigma_array);
    xxfree(min_handle->double_sigma_array);
    
    return(net);
}

static struct fsm *rebuild_machine(struct minimize_handle *min_handle, struct fsm *net) {
  int i,j, group_num, source, target, new_linecount = 0, arccount = 0;
  struct fsm_state *fsm;
  struct p *myp;
  struct e *thise;

  if (net->statecount == min_handle->total_states) {
      return(net);
  }
  fsm = net->states;

  /* We need to make sure state 0 is first in its group */
  /* to get the proper numbering of states */

  if (min_handle->E->group->first_e != min_handle->E) {
    min_handle->E->group->first_e = min_handle->E;
  }

  /* Recycling t_count for group numbering use here */

  group_num = 1;
  myp = min_handle->P;
  while (myp != NULL) {
    myp->count = 0;
    myp = myp->next;
  }

  for (i=0; (fsm+i)->state_no != -1; i++) {
    thise = min_handle->E+((fsm+i)->state_no);
    if (thise->group->first_e == thise) {
      new_linecount++;
      if ((fsm+i)->start_state == 1) {
	thise->group->t_count = 0;
	thise->group->count = 1;
      } else if (thise->group->count == 0) {
	thise->group->t_count = group_num++;
	thise->group->count = 1;
      }
    }
  }

  for (i=0, j=0; (fsm+i)->state_no != -1; i++) {
    thise = min_handle->E+((fsm+i)->state_no);
    if (thise->group->first_e == thise) {
      source = thise->group->t_count;
      target = ((fsm+i)->target == -1) ? -1 : (min_handle->E+((fsm+i)->target))->group->t_count;
      add_fsm_arc(fsm, j, source, (fsm+i)->in, (fsm+i)->out, target, min_handle->finals[(fsm+i)->state_no], (fsm+i)->start_state);
      arccount = ((fsm+i)->target == -1) ? arccount : arccount+1;
      j++;
    }
  }
  
  add_fsm_arc(fsm, j, -1, -1, -1, -1, -1, -1);
  fsm = xxrealloc(fsm,sizeof(struct fsm_state)*(new_linecount+1));
  net->states = fsm;
  net->linecount = j+1;
  net->arccount = arccount;
  net->statecount = min_handle->total_states;
  return(net);
}

static INLINE int refine_states(struct minimize_handle *min_handle, int invstates) {
    int i, selfsplit;
    struct e *thise;
    struct p *tP, *newP = NULL;

  /* 
     1. add inverse(P,a) to table of inverses, disallowing duplicates
     2. first pass on S, touch each state once, increasing P->t_count
     3. for each P where counter != count, split and add to agenda
  */

  /* Inverse to table of inverses */
  selfsplit = 0;

  /* touch and increase P->counter */
  for (i=0; i < invstates; i++) {
    ((min_handle->E+(*(min_handle->temp_move+i)))->group)->t_count++;
    ((min_handle->E+(*(min_handle->temp_move+i)))->group)->inv_t_count += ((min_handle->E+(*(min_handle->temp_move+i)))->inv_count);
    assert((min_handle->E+(*(min_handle->temp_move+i)))->group->t_count <= (min_handle->E+(*(min_handle->temp_move+i)))->group->count);
  }

  /* Split (this is the tricky part) */
  
  for (i=0; i < invstates; i++) {
    
    thise = min_handle->E+*(min_handle->temp_move+i);
    tP = thise->group;
    
    /* Do we need to split?
       if we've touched as many states as there are in the partition
       we don't split */

    if (tP->t_count == tP->count) {
      tP->t_count = 0;
      tP->inv_t_count = 0;
      continue;
    }
    
    if ((tP->t_count != tP->count) && (tP->count > 1) && (tP->t_count > 0)) {      
        
        /* Check if we already split this */
        newP = tP->current_split;
        if (newP == NULL) {
            /* printf("tP [%i] newP [%i]\n",tP->inv_count,tP->inv_t_count); */
            /* Create new group newP */
            min_handle->total_states++;
            if (min_handle->total_states == min_handle->num_states)
                return(1); /* Abort now, machine is already minimal */
            tP->current_split = min_handle->Pnext++;
            newP = tP->current_split;
            newP->first_e = newP->last_e = thise;
            newP->count = 0;
            newP->inv_count = tP->inv_t_count;
            newP->inv_t_count = 0;
            newP->t_count = 0;
            newP->current_split = NULL;
            newP->agenda = NULL;

            /* Add to agenda */
            
            /* If the current block (tP) is on the agenda, we add both back */
            /* to the agenda */
            /* In practice we need only add newP since tP stays where it is */
            /* However, we mark the larger one as not starting the symloop */
            /* from zero */
            if (tP->agenda != NULL) {
                /* Is tP smaller */
                if (tP->inv_count < tP->inv_t_count) {
                    agenda_add(min_handle, newP, 1);
                    tP->agenda->index = 0;
                }
                else {
                    agenda_add(min_handle, newP, 0);
                }
                /* In the event that we're splitting the partition we're currently */
                /* splitting with, we can simply add both new partitions to the agenda */
                /* and break out of the entire sym loop after we're */
                /* done with the current sym and move on with the agenda */
                /* We process the larger one for all symbols */
                /* and the smaller one for only the ones remaining in this symloop */

            } else if (tP == min_handle->current_w) {
                agenda_add(min_handle, ((tP->inv_count < tP->inv_t_count) ? tP : newP),0);
                agenda_add(min_handle, ((tP->inv_count >= tP->inv_t_count) ? tP : newP),1);
                selfsplit = 1;
            } else {
                /* If the block is not on the agenda, we add */
                /* the smaller of tP, newP and start the symloop from 0 */                
                agenda_add(min_handle, (tP->inv_count < tP->inv_t_count ? tP : newP),0);
            }
            /* Add to middle of P-chain */
            newP->next = min_handle->P->next;
            min_handle->P->next = newP;
        }
    
        thise->group = newP;
        newP->count++;
        
        /* need to make tP->last_e point to the last untouched e */
        if (thise == tP->last_e)
            tP->last_e = thise->left;
        if (thise == tP->first_e)
            tP->first_e = thise->right;
        
        /* Adjust links */
        if (thise->left != NULL)
            thise->left->right = thise->right;
        if (thise->right != NULL)
            thise->right->left = thise->left;
        
        if (newP->last_e != thise) {
            newP->last_e->right = thise;
            thise->left = newP->last_e;
            newP->last_e = thise;
        }
    
        thise->right = NULL;
        if (newP->first_e == thise)
            thise->left = NULL;
        
        /* Are we done for this block? Adjust counters */
        if (newP->count == tP->t_count) {
            tP->count = tP->count - newP->count;
            tP->inv_count = tP->inv_count - tP->inv_t_count;
            tP->current_split = NULL;
            tP->t_count = 0;
            tP->inv_t_count = 0;
        }
    }
  }
  /* We return 1 if we just split the partition we were working with */
  return (selfsplit);
}

static void agenda_add(struct minimize_handle *min_handle, struct p *pptr, int start) {

  /* Use FILO strategy here */

  struct agenda *ag;
  //ag = xxmalloc(sizeof(struct agenda));
  ag = min_handle->Agenda_next++;
  if (min_handle->Agenda != NULL)
    ag->next = min_handle->Agenda;
  else
    ag->next = NULL;
  ag->p = pptr;
  ag->index = start;
  min_handle->Agenda = ag;
  pptr->agenda = ag;
}

static void init_PE(struct minimize_handle *min_handle) {
  /* Create two members of P
     (nonfinals,finals)
     and put both of them on the agenda 
  */

  int i;
  struct e *last_f, *last_nonf;
  struct p *nonFP, *FP;
  struct agenda *ag;

  min_handle->mainloop = 1;
  min_handle->memo_table = xxcalloc(min_handle->num_states,sizeof(int));
  min_handle->temp_move = xxcalloc(min_handle->num_states,sizeof(int));
  min_handle->temp_group = xxcalloc(min_handle->num_states,sizeof(int));
  min_handle->Phead = min_handle->P = min_handle->Pnext = xxcalloc(min_handle->num_states+1, sizeof(struct p));
  nonFP = min_handle->Pnext++;
  FP = min_handle->Pnext++;
  nonFP->next = FP;
  nonFP->count = min_handle->num_states-min_handle->num_finals;
  FP->next = NULL;
  FP->count = min_handle->num_finals;
  FP->t_count = 0;
  nonFP->t_count = 0;
  FP->current_split = NULL;
  nonFP->current_split = NULL;
  FP->inv_count = nonFP->inv_count = FP->inv_t_count = nonFP->inv_t_count = 0;
  
  /* How many groups can we put on the agenda? */
  min_handle->Agenda_top = min_handle->Agenda_next = xxcalloc(min_handle->num_states*2, sizeof(struct agenda));
  min_handle->Agenda_head = NULL;

  min_handle->P = NULL;
  min_handle->total_states = 0;

  if (min_handle->num_finals > 0) {
      ag = min_handle->Agenda_next++;
      FP->agenda = ag;
      min_handle->P = FP;
      min_handle->P->next = NULL;
      ag->p = FP;
      min_handle->Agenda_head = ag;
      ag->next = NULL;
      min_handle->total_states++;
  }
  if (min_handle->num_states - min_handle->num_finals > 0) {
      ag = min_handle->Agenda_next++;
      nonFP->agenda = ag;
      ag->p = nonFP;
      ag->next = NULL;
      min_handle->total_states++;
      if (min_handle->Agenda_head != NULL) {
          min_handle->Agenda_head->next = ag;
          min_handle->P->next = nonFP;
          min_handle->P->next->next = NULL;
      } else {
          min_handle->P = nonFP;
          min_handle->P->next = NULL;
          min_handle->Agenda_head = ag;
      }
  }
  
  /* Initialize doubly linked list E */
  min_handle->E = xxcalloc(min_handle->num_states,sizeof(struct e));

  last_f = NULL;
  last_nonf = NULL;
  
  for (i=0; i < min_handle->num_states; i++) {
    if (min_handle->finals[i]) {
      (min_handle->E+i)->group = FP;
      (min_handle->E+i)->left = last_f;
      if (i > 0 && last_f != NULL)
	last_f->right = (min_handle->E+i);
      if (last_f == NULL)
	FP->first_e = (min_handle->E+i);
      last_f = (min_handle->E+i);
      FP->last_e = (min_handle->E+i);
    } else {
      (min_handle->E+i)->group = nonFP;
      (min_handle->E+i)->left = last_nonf;
      if (i > 0 && last_nonf != NULL)
	last_nonf->right = (min_handle->E+i);
      if (last_nonf == NULL)
	nonFP->first_e = (min_handle->E+i);
      last_nonf = (min_handle->E+i);
      nonFP->last_e = (min_handle->E+i);
    }
    (min_handle->E+i)->inv_count = 0;
  }

  if (last_f != NULL)
    last_f->right = NULL;
  if (last_nonf != NULL)
    last_nonf->right = NULL;
}

static int trans_sort_cmp(const void *a, const void *b) {
  return (((const struct trans_list *)a)->inout - ((const struct trans_list *)b)->inout);
}

static void generate_inverse(struct minimize_handle *min_handle, struct fsm *net) {
    
    struct fsm_state *fsm;
    struct trans_array *tptr;
    struct trans_list *listptr;

    int i, source, target, offsetcount, symbol, size;
    fsm = net->states;
    min_handle->trans_array = xxcalloc(net->statecount, sizeof(struct trans_array));
    min_handle->trans_list = xxcalloc(net->arccount, sizeof(struct trans_list));

    /* Figure out the number of transitions each one has */
    for (i=0; (fsm+i)->state_no != -1; i++) {
        if ((fsm+i)->target == -1) {
            continue;
        }
        target = (fsm+i)->target;
        (min_handle->E+target)->inv_count++;
        (min_handle->E+target)->group->inv_count++;
        (min_handle->trans_array+target)->size++;
    }
    offsetcount = 0;
    for (i=0; i < net->statecount; i++) {
        (min_handle->trans_array+i)->transitions = min_handle->trans_list + offsetcount;
        offsetcount += (min_handle->trans_array+i)->size;
    }
    for (i=0; (fsm+i)->state_no != -1; i++) {
        if ((fsm+i)->target == -1) {
            continue;
        }
        symbol = symbol_pair_to_single_symbol(min_handle, (fsm+i)->in,(fsm+i)->out);        
        source = (fsm+i)->state_no;
        target = (fsm+i)->target;
        tptr = min_handle->trans_array + target;
        ((tptr->transitions)+(tptr->tail))->inout = symbol;
        ((tptr->transitions)+(tptr->tail))->source = source;
        tptr->tail++;
    }
    /* Sort arcs */
    for (i=0; i < net->statecount; i++) {
        listptr = (min_handle->trans_array+i)->transitions;
        size = (min_handle->trans_array+i)->size;
        if (size > 1)
            qsort(listptr, size, sizeof(struct trans_list), trans_sort_cmp);
    }
}

static void sigma_to_pairs(struct minimize_handle *min_handle, struct fsm *net) {
    
  int i, j, x, y, z, next_x = 0;
  struct fsm_state *fsm;

  fsm = net->states;
  
  min_handle->epsilon_symbol = -1; 
  min_handle->maxsigma = sigma_max(&net->sigma);

  min_handle->maxsigma++;

  min_handle->single_sigma_array = xxmalloc(2*min_handle->maxsigma*min_handle->maxsigma*sizeof(int));
  min_handle->double_sigma_array = xxmalloc(min_handle->maxsigma*min_handle->maxsigma*sizeof(int));
  
  for (i=0; i < min_handle->maxsigma; i++) {
    for (j=0; j< min_handle->maxsigma; j++) {
      *(min_handle->double_sigma_array+min_handle->maxsigma*i+j) = -1;
    }
  }
  
  /* f(x) -> y,z sigma pair */
  /* f(y,z) -> x simple entry */
  /* if exists f(n) <-> EPSILON, EPSILON, save n */
  /* symbol(x) x>=1 */

  /* Forward mapping: */
  /* *(min_handle->double_sigma_array+min_handle->maxsigma*in+out) */

  /* Backmapping: */
  /* *(min_handle->single_sigma_array+(symbol*2) = in(symbol) */
  /* *(min_handle->single_sigma_array+(symbol*2+1) = out(symbol) */

  /* Table for checking whether a state is final */

  min_handle->finals = xxcalloc(min_handle->num_states, sizeof(_Bool));
  x = 0; min_handle->num_finals = 0;
  net->arity = 1;
  for (i=0; (fsm+i)->state_no != -1; i++) {
    if ((fsm+i)->final_state == 1 && min_handle->finals[(fsm+i)->state_no] != 1) {
      min_handle->num_finals++;
      min_handle->finals[(fsm+i)->state_no] = 1;
    }
    y = (fsm+i)->in;
    z = (fsm+i)->out;
    if (y != z || y == UNKNOWN || z == UNKNOWN)
        net->arity = 2;
    if ((y == -1) || (z == -1))
      continue;
    if (*(min_handle->double_sigma_array+min_handle->maxsigma*y+z) == -1) {
      *(min_handle->double_sigma_array+min_handle->maxsigma*y+z) = x;
      *(min_handle->single_sigma_array+next_x) = y;
      next_x++;
      *(min_handle->single_sigma_array+next_x) = z;
      next_x++;
      if (y == EPSILON && z == EPSILON) {
	min_handle->epsilon_symbol = x;
      }
      x++;
    }
  }
  min_handle->num_symbols = x;
}

static INLINE int symbol_pair_to_single_symbol(struct minimize_handle *min_handle, int in, int out) {
  return(*(min_handle->double_sigma_array+min_handle->maxsigma*in+out));
}

