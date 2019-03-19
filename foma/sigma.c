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

#include <string.h>
#include <stdlib.h>
#include "foma.h"

static int insert_symbol_at(struct sigma *sig, unsigned int i, int num, const char* str)
{
    if (!sig->symbols) {
        sig->symbols = xxcalloc(64, sizeof(struct symbol));
        if (!sig->symbols)
            return -1;
        sig->size = 0;
        sig->memsize = 64;
    } else if (sig->memsize <= sig->size) {
        sig->memsize = (unsigned int)(sig->memsize * 2);
        sig->symbols = xxrealloc(sig->symbols, sig->memsize * sizeof(struct symbol));
        if (!sig->symbols)
            return -1;
    }

    if (sig->size <= i) {
        sig->symbols[sig->size].number = num;
        sig->symbols[sig->size].symbol = xxstrdup(str);
    } else {
        memmove(sig->symbols + i + 1, sig->symbols + i, (sig->size - i) * sizeof(struct symbol));
        sig->symbols[i].number = num;
        sig->symbols[i].symbol = xxstrdup(str);
    }
    sig->size++;

    return num;
}

void remove_symbol_from(struct sigma *sig, unsigned int i)
{
    xxfree(sig->symbols[i].symbol);
    memmove(sig->symbols + i, sig->symbols + i + 1, sizeof(struct symbol) * (sig->size - (i + 1)));
    sig->size--;
}

static unsigned int begin_index_of_non_special_symbol(struct sigma *sig)
{
    for (unsigned int i = 0; i < sig->size; ++i) {
        if (IDENTITY < sig->symbols[i].number)
            return i;
    }

    return 0;
}

int sigma_remove(const char *symbol, struct sigma *sigma) {
  if (!sigma || !symbol)
      return -1;

  for (unsigned int i = 0; i < sigma->size; ++i) {
      if (strcmp(sigma->symbols[i].symbol, symbol) == 0) {
          remove_symbol_from(sigma, i);
          return i;
      }
  }
  return -1;
}

int sigma_remove_num(int num, struct sigma *sigma) {
  if (!sigma || num == -1)
      return -1;

  for (unsigned int i = 0; i < sigma->size; ++i) {
      if (sigma->symbols[i].number == num) {
          remove_symbol_from(sigma, i);
          return i;
      }
  }
  return -1;
}

int sigma_add_special (int symbol, struct sigma *sigma) {
    unsigned int i;

    if (!sigma)
        return -1;

    /* Insert special symbols pre-sorted */
    for (i = 0; i < sigma->size; ++i) {
        if (symbol <= sigma->symbols[i].number)
            break;
    }

    if (symbol == EPSILON)
        return insert_symbol_at(sigma, i, EPSILON, "@_EPSILON_SYMBOL_@");
    else if (symbol == IDENTITY)
        return insert_symbol_at(sigma, i, IDENTITY, "@_IDENTITY_SYMBOL_@");
    else if (symbol == UNKNOWN)
        return insert_symbol_at(sigma, i, UNKNOWN, "@_UNKNOWN_SYMBOL_@");
    else
        return -1;
}

/* WARNING: this function will indeed add a symbol to sigma */
/* but it's up to the user to sort the sigma (affecting arc numbers in the network) */
/* before merge_sigma() is ever called */

int sigma_add (const char *symbol, struct sigma *sigma) {
    int num;

    if (!sigma || !symbol)
        return -1;

    /* Insert special symbols pre-sorted */
    if (strcmp(symbol, "@_EPSILON_SYMBOL_@") == 0)
        return sigma_add_special(EPSILON, sigma);
    else if (strcmp(symbol,"@_IDENTITY_SYMBOL_@") == 0)
        return sigma_add_special(IDENTITY, sigma);
    else if (strcmp(symbol,"@_UNKNOWN_SYMBOL_@") == 0)
        return sigma_add_special(UNKNOWN, sigma);

    /* Insert non-special in any order */
    if (sigma->size == 0)
        return insert_symbol_at(sigma, 0, 3, symbol);

    num = sigma->symbols[sigma->size - 1].number + 1;
    if (num <= IDENTITY)
        num = 3;

    if (insert_symbol_at(sigma, sigma->size, num, symbol) < 0)
        return -1;

    return num;
}

/* Remove symbols that are never used from sigma and renumber   */
/* The variable force controls whether to remove even though    */
/* @ or ? is present                                            */
/* If force == 1, unused symbols are always removed regardless  */

void sigma_cleanup (struct fsm *net, int force) {
    int i,j,maxsigma,*attested;
    struct fsm_state *fsm;
    struct symbol *syms = net->sigma.symbols;

    if (force == 0) {
        if (sigma_find_number(IDENTITY, &net->sigma) != -1)
            return;
        if (sigma_find_number(UNKNOWN, &net->sigma) != -1)
            return;
    }

    maxsigma = sigma_max(&net->sigma);
    if (maxsigma < 0) { return; }
    attested = xxmalloc(sizeof(int)*(maxsigma+1));
    memset(attested, 0, sizeof(int) * (maxsigma + 1));

    fsm = net->states;
    for (i=0; (fsm+i)->state_no != -1; i++) {
        if ((fsm+i)->in >=0)
            *(attested+(fsm+i)->in) = 1;
        if ((fsm+i)->out >=0)
            *(attested+(fsm+i)->out) = 1;
    }
    for (i = 3, j = 3; i <= maxsigma; ++i) {
        if (*(attested+i)) {
            *(attested+i) = j;
            j++;
        }
    }
    for (i=0; (fsm+i)->state_no != -1; i++) {
        if ((fsm+i)->in > 2)
            (fsm+i)->in = *(attested+(fsm+i)->in);
        if ((fsm+i)->out > 2)
            (fsm+i)->out = *(attested+(fsm+i)->out);
    }

    for (i = 0; i < net->sigma.size; ++i) {
        if (!attested[syms[i].number])
            remove_symbol_from(&net->sigma, i);
        else if (syms[i].number > 2)
            syms[i].number = attested[syms[i].number];
    }
    xxfree(attested);
    return;
}

int sigma_max(struct sigma *sigma) {
  int max = -1;
  struct symbol *syms;

  if (sigma == NULL)
      return -1;

  syms = sigma->symbols;
  for (unsigned int i = 0; i < sigma->size; ++i) {
      if (max < syms[i].number)
          max = syms[i].number;
  }

  return max;
}

struct fsm_sigma_list *sigma_to_list(struct sigma *sigma) {
    struct fsm_sigma_list *sl;
    sl = xxcalloc(sigma_max(sigma)+1,sizeof(struct fsm_sigma_list));
    for (unsigned int i = 0; i < sigma->size; ++i)
        sl[sigma->symbols[i].number].symbol = sigma->symbols[i].symbol;
    return sl;
}

int sigma_add_number(struct sigma *sigma, const char *symbol, int number) {
    if (!sigma)
        return -1;

    return insert_symbol_at(sigma, sigma->size, number, symbol);
}

int sigma_find_number(int number, const struct sigma *sigma) {
    if (sigma == NULL)
        return -1;

    for (unsigned int i = 0; i < sigma->size; ++i) {
        if (sigma->symbols[i].number == number)
            return number;
    }

    return -1;
}
char *sigma_string(int number, struct sigma *sigma) {
    if (sigma == NULL)
        return NULL;

    for (unsigned int i = 0; i < sigma->size; ++i) {
        if (sigma->symbols[i].number == number)
            return sigma->symbols[i].symbol;
    }

    return NULL;
}

/* Substitutes string symbol for sub in sigma */
/* no check for duplicates                    */
int sigma_substitute(const char *symbol, const char *sub, struct sigma *sigma) {
    if (!sigma || !symbol || !sub)
        return -1;

    for (unsigned int i = 0; i < sigma->size; ++i) {
        if (strcmp(sigma->symbols[i].symbol, symbol) == 0) {
            xxfree(sigma->symbols[i].symbol);
            sigma->symbols[i].symbol = strdup(sub);
            return sigma->symbols[i].number;
        }
    }

    return -1;
}

int sigma_find(const char *symbol, const struct sigma *sigma) {
    if (!sigma || !symbol)
        return -1;

    for (unsigned int i = 0; i < sigma->size; ++i) {
        if (strcmp(sigma->symbols[i].symbol, symbol) == 0)
            return sigma->symbols[i].number;
    }

    return -1;
}

int symbol_comp(const struct symbol *a, const struct symbol *b) {
  return(strcmp(a->symbol, b->symbol));
}

struct sigma sigma_copy(const struct sigma *sigma) {
    struct sigma copied = sigma_create();

    if (sigma == NULL)
        return copied;

    copied.size = sigma->size;
    copied.memsize = sigma->memsize;
    copied.symbols = xxcalloc(copied.memsize, sizeof(struct symbol));

    for (unsigned int i = 0; i < copied.size; ++i) {
        copied.symbols[i].number = sigma->symbols[i].number;
        if (sigma->symbols[i].symbol)
            copied.symbols[i].symbol = xxstrdup(sigma->symbols[i].symbol);
        else
            copied.symbols[i].symbol = NULL;
    }

    return copied;
}

/* Assigns a consecutive numbering to symbols in sigma > IDENTITY */
/* and sorts the sigma based on the symbol string contents        */

int sigma_sort(struct fsm *net) {
  int(*comp)() = symbol_comp;
  int max, *replacearray;
  struct fsm_state *fsm_state = net->states;
  struct symbol *syms = net->sigma.symbols;
  int start = begin_index_of_non_special_symbol(&net->sigma);

  if (net->sigma.size == 0)
      return 1;

  qsort(syms + start, net->sigma.size - start, sizeof(struct symbol), comp);

  max = sigma_max(&net->sigma);
  replacearray = xxmalloc(sizeof(int) * (max + 3));

  if (0 < net->sigma.size) {
      unsigned int i = start;
      int new_num = 3;
      while (i < net->sigma.size - 1) {
          replacearray[syms[i].number] = new_num;
          if (strcmp(syms[i].symbol, syms[i + 1].symbol) == 0) {
              remove_symbol_from(&net->sigma, i);
          } else {
              syms[i].number = new_num;
              i++;
              new_num++;
          }
      }

      if (i < net->sigma.size)
          replacearray[syms[i].number] = new_num;
  }

  /* Replace arcs */
  for(unsigned int i = 0; fsm_state[i].state_no != -1; ++i) {
      if (fsm_state[i].in > IDENTITY)
          fsm_state[i].in = replacearray[fsm_state[i].in];
      if (fsm_state[i].out > IDENTITY)
          fsm_state[i].out = replacearray[fsm_state[i].out];
  }

  xxfree(replacearray);

  return 1;
}

struct sigma sigma_create() {
    struct sigma sig = {.symbols = NULL, .size = 0, .memsize = 0};
    return sig;
}


struct merge_symbol {
    int number;
    const char *symbol;
    unsigned char presence; /* 1 = in net 1, 2 = in net 2, 3 = in both */
};

struct merge_sigma {
  struct merge_symbol *symbols;
  unsigned int size;
  unsigned int memsize;
};

static void append_to_mergesigma(struct merge_sigma *msigma, struct symbol sym, short presence) {
    int number = 0;

    if (msigma->size == 0) {
        number = 2;
    } else {
        number = msigma->symbols[msigma->size - 1].number;
    }

    if (sym.number < 3) {
        msigma->symbols[msigma->size].number = sym.number;
    } else {
        if (number < 3)
            number = 2;
        msigma->symbols[msigma->size].number = number + 1;
    }
    msigma->symbols[msigma->size].symbol = sym.symbol;
    msigma->symbols[msigma->size].presence = presence;
    msigma->size++;
}

static struct sigma copy_mergesigma(struct merge_sigma *msigma) {
    struct sigma sigma = { .symbols = NULL, .size = msigma->size, .memsize = msigma->memsize };
    struct merge_symbol *msyms = msigma->symbols;

    sigma.symbols = xxcalloc(sigma.memsize, sizeof(struct symbol));

    for (unsigned int i = 0; i < msigma->size; ++i) {
        sigma.symbols[i].number = msyms[i].number;
        sigma.symbols[i].symbol = xxstrdup(msyms[i].symbol);
    }

    return sigma;
}

static struct fsm_state *create_merged_states(struct fsm_state *old_states, struct merge_sigma msigma, int presence)
{
    int net_unk = 0;
    int net_adds = 0;
    int net_lines = find_arccount(old_states);
    struct fsm_state *new_states;
    int j = 0;

    for (unsigned int i = 0; i < msigma.size; ++i) {
        if (msigma.symbols[i].presence == presence)
            net_unk++;
    }

    for (unsigned int i = 0; old_states[i].state_no != -1; i++) {
        if ((old_states+i)->in == IDENTITY)
            net_adds += net_unk;
        if (((old_states+i)->in == UNKNOWN) && ((old_states+i)->out != UNKNOWN))
            net_adds += net_unk;
        if (((old_states+i)->out == UNKNOWN) && ((old_states+i)->in != UNKNOWN))
            net_adds += net_unk;
        if (((old_states+i)->in == UNKNOWN) && ((old_states+i)->out == UNKNOWN))
            net_adds += net_unk*net_unk - net_unk + 2*net_unk;
    }

    new_states = xxmalloc(sizeof(struct fsm_state)*(net_adds+net_lines+1));
    j = 0;
    for (int i = 0; (old_states+i)->state_no != -1; i++) {

        if ((old_states+i)->in == IDENTITY) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
            for (unsigned int m = 0; m < msigma.size; ++m) {
                if ((msigma.symbols[m].presence == presence) && (msigma.symbols[m].number > IDENTITY)) {
                    add_fsm_arc(new_states, j, (old_states+i)->state_no,
                        msigma.symbols[m].number, msigma.symbols[m].number, (old_states+i)->target,
                        (old_states+i)->final_state, (old_states+i)->start_state);
                    j++;
                }
            }
        }

        else if ((old_states+i)->in == UNKNOWN && (old_states+i)->out != UNKNOWN) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
            for (unsigned int m = 0; m < msigma.size; ++m) {
                if ((msigma.symbols[m].presence == presence) && (msigma.symbols[m].number > IDENTITY)) {
                    add_fsm_arc(new_states, j, (old_states+i)->state_no,
                        msigma.symbols[m].number, (old_states+i)->out, (old_states+i)->target,
                        (old_states+i)->final_state, (old_states+i)->start_state);
                    j++;
                }
            }
        }

        else if (((old_states+i)->in != UNKNOWN) && ((old_states+i)->out == UNKNOWN)) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
            for (unsigned int m = 0; m < msigma.size; ++m) {
                if ((msigma.symbols[m].presence == presence) && (msigma.symbols[m].number > IDENTITY)) {
                    add_fsm_arc(new_states, j, (old_states+i)->state_no,
                        (old_states+i)->in, msigma.symbols[m].number, (old_states+i)->target,
                        (old_states+i)->final_state, (old_states+i)->start_state);
                    j++;
                }
            }
        }

        /* Replace ?:? with ?:[all unknowns] [all unknowns]:? and [all unknowns]:[all unknowns] where a != b */
        else if ((old_states+i)->in == UNKNOWN && (old_states+i)->out == UNKNOWN) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
            for (unsigned int m1 = 0; m1 < msigma.size; ++m1) {
                for (unsigned int m2 = 0; m2 < msigma.size; ++m2) {
                    if (msigma.symbols[m1].number != msigma.symbols[m2].number
                        && ((msigma.symbols[m1].presence == presence && msigma.symbols[m2].presence == presence && msigma.symbols[m1].number > IDENTITY && msigma.symbols[m2].number > IDENTITY)
                            || (msigma.symbols[m1].number == UNKNOWN && msigma.symbols[m2].number > IDENTITY && msigma.symbols[m2].presence == presence)
                            || (msigma.symbols[m2].number == UNKNOWN && msigma.symbols[m1].number > IDENTITY && msigma.symbols[m1].presence == presence))) {
                        add_fsm_arc(new_states, j, (old_states+i)->state_no,
                            msigma.symbols[m1].number, msigma.symbols[m2].number, (old_states+i)->target,
                            (old_states+i)->final_state, (old_states+i)->start_state);
                        j++;
                    }
                }
            }
        }

        /* Simply copy arcs that are not IDENTITY or UNKNOWN */
        else if (((old_states+i)->in > IDENTITY || (old_states+i)->in == EPSILON)
            && ((old_states+i)->out > IDENTITY || (old_states+i)->out == EPSILON)) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
        }

        else if ((old_states+i)->in == -1) {
            add_fsm_arc(new_states, j, (old_states+i)->state_no,
                (old_states+i)->in, (old_states+i)->out, (old_states+i)->target,
                (old_states+i)->final_state, (old_states+i)->start_state);
            j++;
        }
    }

    add_fsm_arc(new_states, j, -1, -1, -1, -1, -1, -1);

    return new_states;
}

void fsm_merge_sigma(struct fsm *net1, struct fsm *net2) {

  struct sigma new_sigma_1, new_sigma_2;
  struct merge_sigma msigma;
  struct fsm_state *fsm_state;
  int i, j, end_1 = 0, end_2 = 0, sigmasizes, *mapping_1, *mapping_2, equal = 1, unknown_1 = 0, unknown_2 = 0;
  struct symbol *syms1, *syms2;

  i = sigma_find(".#.", &net1->sigma);
  j = sigma_find(".#.", &net2->sigma);
  if (i != -1 && j == -1) {
      sigma_add(".#.", &net2->sigma);
      sigma_sort(net2);
  }
  if (j != -1 && i == -1) {
      sigma_add(".#.", &net1->sigma);
      sigma_sort(net1);
  }

  sigmasizes = net1->sigma.size + net2->sigma.size;

  mapping_1 = xxmalloc(sizeof(int)*(sigmasizes+3));
  mapping_2 = xxmalloc(sizeof(int)*(sigmasizes+3));

  /* Fill mergesigma */

  msigma.size = 0;
  msigma.memsize = sigmasizes;
  msigma.symbols = xxcalloc(msigma.memsize, sizeof(struct merge_symbol));
  memset(msigma.symbols, 0, sizeof(struct merge_symbol) * msigma.memsize);

  syms1 = net1->sigma.symbols;
  syms2 = net2->sigma.symbols;

  /* Loop over sigma 1, sigma 2 */
  for (unsigned int s1 = 0, s2 = 0, m = 0;;) {
      if (net1->sigma.size <= s1)
          end_1 = 1;
      if (net2->sigma.size <= s2)
          end_2 = 1;
      if (end_1 && end_2)
          break;

      if (end_2) {
          append_to_mergesigma(&msigma, syms1[s1], 1);
          mapping_1[syms1[s1].number] = msigma.symbols[m].number;
          m++;
          s1++;
          equal = 0;
      } else if (end_1) {
          append_to_mergesigma(&msigma, syms2[s2], 2);
          mapping_2[syms2[s2].number] = msigma.symbols[m].number;
          m++;
          s2++;
          equal = 0;
      }

      /* Both alive */

      /* 1 or 2 contains special characters */

      else if ((syms1[s1].number <= IDENTITY) || (syms2[s2].number <= IDENTITY)) {
          /* Treating zeros or unknowns */

          unknown_1 |= (syms1[s1].number == UNKNOWN) || (syms1[s1].number == IDENTITY);
          unknown_2 |= (syms2[s2].number == UNKNOWN) || (syms2[s2].number == IDENTITY);

          if (syms1[s1].number == syms2[s2].number) {
              append_to_mergesigma(&msigma, syms1[s1], 3);
              m++;
              s1++;
              s2++;
          } else if (syms1[s1].number < syms2[s2].number) {
              append_to_mergesigma(&msigma, syms1[s1], 1);
              m++;
              s1++;
              equal = 0;
          } else {
              append_to_mergesigma(&msigma, syms2[s2], 2);
              m++;
              s2++;
              equal = 0;
          }
      } else if (strcmp(syms1[s1].symbol, syms2[s2].symbol) == 0) {  /* Both contain non-special chars */
          append_to_mergesigma(&msigma, syms1[s1], 3);
          /* Add symbol numbers to mapping */
          mapping_1[syms1[s1].number] = msigma.symbols[m].number;
          mapping_2[syms2[s2].number] = msigma.symbols[m].number;

          m++;
          s1++;
          s2++;
      } else if (strcmp(syms1[s1].symbol, syms2[s2].symbol) < 0) {
          append_to_mergesigma(&msigma, syms1[s1], 1);
          mapping_1[syms1[s1].number] = msigma.symbols[m].number;
          m++;
          s1++;
          equal = 0;
      } else {
          append_to_mergesigma(&msigma, syms2[s2], 2);
          mapping_2[syms2[s2].number] = msigma.symbols[m].number;
          m++;
          s2++;
          equal = 0;
      }
  }

  /* Go over both net1 and net2 and replace arc numbers with new mappings */

  fsm_state = net1->states;
  for (i=0; (fsm_state+i)->state_no != -1; i++) {
    if ((fsm_state+i)->in > 2)
      (fsm_state+i)->in = *(mapping_1+(fsm_state+i)->in);
    if ((fsm_state+i)->out > 2)
      (fsm_state+i)->out = *(mapping_1+(fsm_state+i)->out);
  }
  fsm_state = net2->states;
  for (i=0; (fsm_state+i)->state_no != -1; i++) {
    if ((fsm_state+i)->in > 2)
      (fsm_state+i)->in = *(mapping_2+(fsm_state+i)->in);
    if ((fsm_state+i)->out > 2)
      (fsm_state+i)->out = *(mapping_2+(fsm_state+i)->out);
  }

  /* Expand on ?, ?:x, y:? */

  if (unknown_1 && !equal) {
      /* Expand net 1 */
      struct fsm_state *new_states = create_merged_states(net1->states, msigma, 2);
      xxfree(net1->states);
      net1->states = new_states;
  }

  if (unknown_2 && !equal) {
      /* Expand net 2 */
      struct fsm_state *new_states = create_merged_states(net2->states, msigma, 1);
      xxfree(net2->states);
      net2->states = new_states;
  }

  /* Copy mergesigma to net1, net2 */

  new_sigma_1 = copy_mergesigma(&msigma);
  new_sigma_2 = copy_mergesigma(&msigma);

  fsm_sigma_destroy(&net1->sigma);
  fsm_sigma_destroy(&net2->sigma);

  net1->sigma = new_sigma_1;
  net2->sigma = new_sigma_2;

  xxfree(mapping_1);
  xxfree(mapping_2);
  xxfree(msigma.symbols);
}

