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
#include <string.h>
#include "foma/foma_r.h"

// Lower(X) puts X on output tape (may also be represented by @ID@ on input tape)
// Upper(X) puts X on input tape
// Unrewritten(X) X on input tape, not rewritten (aligned with @O@ symbols)
// NotContain(X) MT does not contain MT configuration X

// Boundary: every MT word begins and ends with boundary, i.e. the @#@ symbol on the input tape, output tape, and relevant semantic symbols

/*

       [ @O@  ]  [ @I[@        ] [ @I@         ] [ @I]@        ] [ @I[]@       ]
       [ @0@  ]  [ @#0001@     ] [ @#0001@     ] [ @#0001@     ] [ @#0001@     ]
       [ @#@  ]  [ ANY|@0@     ] [ ANY|@0@     ] [ ANY|@0@     ] [ ANY|@0@     ]
       [ @ID@ ]  [ ANY|@ID|@0@ ] [ ANY|@ID|@0@ ] [ ANY|@ID|@0@ ] [ ANY|@ID|@0@ ]


*/
/* Special symbols used:
   @0@    Epsilon
   @O@    Outside rewrite
   @I@    Inside rewrite
   @I[@   Beginning of rewrite
   @I[]@  Beginning and end of rewrite
   @I]@   End of rewrite
   @ID@   Identity symbol (= repeat symbol on previous tape at this position)
   @#@    Boundary (the symbol .#. is mapped to this in contexts before the compilation procedure)
   @#X@   X = rule number (one for each rule, starting with @#0001@)
*/

struct rewrite_batch {

    struct rewrite_set *rewrite_set;
    struct fsm *Rulenames;
    struct fsm *ISyms;
    struct fsm *ANY;
    struct fsm *IOpen;
    struct fsm *IClose;
    struct fsm *ITape;
    struct fsm *Any4Tape;
    struct fsm *Epextend;
    int num_rules;
    char (*namestrings)[8];

};

char *specialsymbols[] = {"@0@","@O@","@I@","@I[@","@I[]@","@I]@","@ID@","@#@", NULL};

void rewrite_add_special_syms(struct rewrite_batch *rb, struct fsm *net);
struct fsm *rewrite_upper(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper);
struct fsm *rewrite_lower(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lower);
struct fsm *rewrite_two_level(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rightside);
struct fsm *rewr_context_restrict(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *X, struct fsmcontexts *LR);
struct fsm *rewr_contains(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang);
struct fsm *rewr_unrewritten(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang);
struct fsm *rewr_notleftmost(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number, int arrow_type);
struct fsm *rewr_notshortest(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number);
struct fsm *rewr_notlongest(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number, int arrow_type);
struct fsm *rewrite_tape_m_to_n_of_k(struct build_handle *b_handle, struct fsm *lang, int m, int n, int k);
struct fsm *rewrite_cp(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper, struct fsm *lower, int rule_number);
struct fsm *rewrite_cp_transducer(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *t, int rule_number);
struct fsm *rewrite_cp_markup(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper, struct fsm *lower1, struct fsm *lower2, int rule_number);
struct fsm *rewrite_epextend(struct build_handle *b_handle, struct rewrite_batch *rb);
struct fsm *rewrite_any_4tape(struct build_handle *b_handle, struct rewrite_batch *rb);
struct fsm *rewrite_itape(struct rewrite_batch *rb);
void rewrite_cleanup(struct rewrite_batch *rb);


struct fsm *fsm_rewrite(struct build_handle *b_handle, struct rewrite_set *all_rules) {
    struct rewrite_batch *rb;
    struct rewrite_set *ruleset;
    struct fsmrules *rules;
    struct fsmcontexts *contexts;
    struct fsm *RuleCP, *Base, *Boundary, *Outside, *CP, *C, *LeftExtend, *RightExtend, *Center;
    int i, num_rules, rule_number, dir;
    /* Count parallel rules */
    for (ruleset = all_rules, num_rules = 0; ruleset != NULL; ruleset = ruleset->next) {
	for (rules = ruleset->rewrite_rules; rules != NULL; rules = rules->next) {
	     num_rules++;
	 }
    }

    rb = xxcalloc(1, sizeof(struct rewrite_batch));
    rb->rewrite_set = all_rules;
    rb->num_rules = num_rules;
    rb->namestrings = xxmalloc(sizeof *rb->namestrings * num_rules);
    for (i = 0; i < rb->num_rules; i++) {
	sprintf(rb->namestrings[i], "@#%04i@", i+1);
    }

    rb->ISyms = fsm_minimize(b_handle, fsm_union(fsm_symbol("@I@"), fsm_union(fsm_symbol("@I[]@"), fsm_union(fsm_symbol("@I[@"), fsm_symbol("@I]@")))));
    rb->Rulenames = fsm_empty_set();
    for (i = 1; i <= num_rules; i++) {
	rb->Rulenames = fsm_minimize(b_handle, fsm_union(rb->Rulenames, fsm_symbol(rb->namestrings[i-1])));
    }
    rb->ANY = fsm_identity();
    rewrite_add_special_syms(rb, rb->ANY);

    /* Add auxiliary symbols to all alphabets */
    for (ruleset = all_rules; ruleset != NULL; ruleset = ruleset->next) {
        for (rules = ruleset->rewrite_rules; rules != NULL; rules = rules->next) {
            rewrite_add_special_syms(rb, rules->left);
	    rewrite_add_special_syms(rb, rules->right);
            rewrite_add_special_syms(rb, rules->right2);
        }
        for (contexts = ruleset->rewrite_contexts; contexts != NULL; contexts = contexts->next) {
            rewrite_add_special_syms(rb, contexts->left);
            rewrite_add_special_syms(rb, contexts->right);
        }
    }
    /* Get cross-product of every rule, according to its type */
    RuleCP = fsm_empty_set();
    for (ruleset = all_rules, rule_number = 1; ruleset != NULL; ruleset = ruleset->next) {
	dir = ruleset->rule_direction;
        for (rules = ruleset->rewrite_rules; rules != NULL; rules = rules->next) {
	    if (rules->right == NULL) {
		/* T(x)-type rule */
		CP = rewrite_cp_transducer(b_handle, rb, fsm_copy(rules->left), rule_number);
		rules->cross_product = fsm_copy(CP);
		rules->right = fsm_minimize(b_handle, fsm_lower(b_handle, fsm_copy(rules->left)));
		rules->left = fsm_minimize(b_handle, fsm_upper(b_handle, fsm_copy(rules->left)));
		rewrite_add_special_syms(rb, rules->right);
		rewrite_add_special_syms(rb, rules->left);
	    } else if (rules->right2 == NULL) {
		/* Regular rewrite rule */
		CP = rewrite_cp(b_handle, rb, fsm_copy(rules->left), fsm_copy(rules->right), rule_number);
		rules->cross_product = fsm_copy(CP);
	    } else {
		/* A -> B ... C -type rule */
		CP = rewrite_cp_markup(b_handle, rb, fsm_copy(rules->left), fsm_copy(rules->right), fsm_copy(rules->right2), rule_number);
		rules->cross_product = fsm_copy(CP);
	    }
	    RuleCP = fsm_minimize(b_handle, fsm_union(RuleCP, CP));
	    rule_number++;
	}
    }

    /* Create Base language */
    Boundary = fsm_parse_regex("\"@O@\" \"@0@\" \"@#@\" \"@ID@\"", NULL, NULL);
    Outside = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_symbol("@O@"), fsm_concat(b_handle, fsm_symbol("@0@"), fsm_concat(b_handle, fsm_copy(rb->ANY), fsm_symbol("@ID@")))));
    Base = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_copy(Boundary), fsm_concat(b_handle, fsm_kleene_star(b_handle, fsm_union(RuleCP, Outside)), fsm_copy(Boundary))));
    fsm_destroy(Boundary);
    for (ruleset = all_rules, rule_number = 1; ruleset != NULL; ruleset = ruleset->next) {
	dir = ruleset->rule_direction;
	/* Replace all context spec with Upper/Lower, depending on rule_direction */
	for (contexts = ruleset->rewrite_contexts; contexts != NULL; contexts = contexts->next) {
	    switch(dir) {
	    case OP_UPWARD_REPLACE:
		contexts->cpleft = rewrite_upper(b_handle, rb, fsm_copy(contexts->left));
		contexts->cpright = rewrite_upper(b_handle, rb, fsm_copy(contexts->right));
		break;
	    case OP_RIGHTWARD_REPLACE:
		contexts->cpleft = rewrite_lower(b_handle, rb, fsm_copy(contexts->left));
		contexts->cpright = rewrite_upper(b_handle, rb, fsm_copy(contexts->right));
		break;
	    case OP_LEFTWARD_REPLACE:
		contexts->cpleft = rewrite_upper(b_handle, rb, fsm_copy(contexts->left));
		contexts->cpright = rewrite_lower(b_handle, rb, fsm_copy(contexts->right));
		break;
	    case OP_DOWNWARD_REPLACE:
		contexts->cpleft = rewrite_lower(b_handle, rb, fsm_copy(contexts->left));
		contexts->cpright = rewrite_lower(b_handle, rb, fsm_copy(contexts->right));
		break;
	    case OP_TWO_LEVEL_REPLACE:
		contexts->cpleft = rewrite_two_level(b_handle, rb, fsm_copy(contexts->left), 0);
		contexts->cpright = rewrite_two_level(b_handle, rb, fsm_copy(contexts->right), 1);
		break;
	    }
	}
        for (rules = ruleset->rewrite_rules; rules != NULL; rules = rules->next) {
	    /* Just the rule center w/ number without CP() contests */
	    /* Actually, maybe better to include CP(U,L) in this, very slow with e.g. a -> a || _ b^15 */
	    if (rules->arrow_type & ARROW_DOTTED) {
		/* define EP Tape1of4("@O@") | [ Tape1of4("@I[@" "@I@"* "@I]@" | "@I[]@") & Tape3of4(~["@0@"*]) ] ; */
		/* Additional constraint: 0->x is only allowed between EP _ EP */
		/* The left and right sides can be checked separately */
		/* ~[?* Center ~[EP ?*]] & ~[~[?* EP] Center ?*] */
		Center = fsm_copy(rules->cross_product);
		Base = fsm_intersect(b_handle, fsm_intersect(b_handle, Base, fsm_complement(b_handle, fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), fsm_concat(b_handle, fsm_copy(Center), fsm_complement(b_handle, fsm_concat(b_handle, rewrite_epextend(b_handle, rb), rewrite_any_4tape(b_handle, rb))))))), fsm_complement(b_handle, fsm_concat(b_handle, fsm_complement(b_handle, fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), rewrite_epextend(b_handle, rb))), fsm_concat(b_handle, fsm_copy(Center), rewrite_any_4tape(b_handle, rb)))));
	    }
	    if (ruleset->rewrite_contexts) {
		Base = fsm_intersect(b_handle, Base, rewr_context_restrict(b_handle, rb, rules->cross_product, ruleset->rewrite_contexts));
	    }
	    /* Determine C (based on rule type) */
	    C = fsm_empty_set();
	    if ((rules->arrow_type & ARROW_RIGHT) && !(rules->arrow_type & ARROW_OPTIONAL)) {
		C = fsm_union(C, rewr_unrewritten(b_handle, rb, fsm_minimize(b_handle, fsm_minus(b_handle, fsm_copy(rules->left), fsm_empty_string()))));
	    }
	    if ((rules->arrow_type & ARROW_LEFT) && !(rules->arrow_type & ARROW_OPTIONAL)) {
		C = fsm_union(C, rewr_unrewritten(b_handle, rb, fsm_minimize(b_handle, fsm_minus(b_handle, fsm_copy(rules->right), fsm_empty_string()))));
	    }
	    if (rules->arrow_type & ARROW_LONGEST_MATCH) {
		if (rules->arrow_type & ARROW_RIGHT) {
		    C = fsm_union(C, rewr_notleftmost(b_handle, rb, rewrite_upper(b_handle, rb, fsm_copy(rules->left)), rule_number, rules->arrow_type));
		    C = fsm_union(C, rewr_notlongest(b_handle, rb, rewrite_upper(b_handle, rb, fsm_copy(rules->left)), rule_number, rules->arrow_type));
		}
		if (rules->arrow_type & ARROW_LEFT) {
		    C = fsm_union(C, rewr_notleftmost(b_handle, rb, rewrite_lower(b_handle, rb, fsm_copy(rules->right)), rule_number, rules->arrow_type));
		    C = fsm_union(C, rewr_notlongest(b_handle, rb, rewrite_lower(b_handle, rb, fsm_copy(rules->right)), rule_number, rules->arrow_type));
		}
	    }
	    if (rules->arrow_type & ARROW_SHORTEST_MATCH) {
		if (rules->arrow_type & ARROW_RIGHT) {		
		    C = fsm_union(C, rewr_notleftmost(b_handle, rb, rewrite_upper(b_handle, rb, fsm_copy(rules->left)), rule_number, rules->arrow_type));
		    C = fsm_union(C, rewr_notshortest(b_handle, rb, rewrite_upper(b_handle, rb, fsm_copy(rules->left)), rule_number));
		}
		if (rules->arrow_type & ARROW_LEFT) {
		    C = fsm_union(C, rewr_notleftmost(b_handle, rb, rewrite_lower(b_handle, rb, fsm_copy(rules->right)), rule_number, rules->arrow_type));
		    C = fsm_union(C, rewr_notshortest(b_handle, rb, rewrite_lower(b_handle, rb, fsm_copy(rules->right)), rule_number));
		}
	    }
	    if (!ruleset->rewrite_contexts) {
		if (rules->arrow_type & ARROW_DOTTED && !(rules->arrow_type & ARROW_OPTIONAL)) {
		    Base = fsm_minus(b_handle, Base, rewr_contains(b_handle, rb, fsm_concat(b_handle, rewrite_epextend(b_handle, rb), rewrite_epextend(b_handle, rb))));
		} else {
		    Base = fsm_minus(b_handle, Base, rewr_contains(b_handle, rb, fsm_copy(C)));
		}
	    }
	    for (contexts = ruleset->rewrite_contexts; contexts != NULL; contexts = contexts->next) {
		/* Constraints: running intersect w/ Base */
		/* NotContain(LC [Unrewritten|LM|...] RC) */
		if (rules->arrow_type & ARROW_DOTTED && !(rules->arrow_type & ARROW_OPTIONAL)) {
		    /* Extend left and right */
		    LeftExtend = fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), fsm_copy(contexts->cpleft)), fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), rewrite_epextend(b_handle, rb))));
		    RightExtend = fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_concat(b_handle, rewrite_epextend(b_handle, rb), rewrite_any_4tape(b_handle, rb)), fsm_concat(b_handle, fsm_copy(contexts->cpright), rewrite_any_4tape(b_handle, rb))));
		    Base = fsm_minus(b_handle, Base, rewr_contains(b_handle, rb, fsm_minimize(b_handle, fsm_concat(b_handle, LeftExtend, RightExtend))));
		} else {
		    Base = fsm_minus(b_handle, Base, rewr_contains(b_handle, rb, fsm_concat(b_handle, fsm_copy(contexts->cpleft), fsm_concat(b_handle, fsm_copy(C), fsm_copy(contexts->cpright)))));
		}
	    }
	    rule_number++;
	    fsm_destroy(C);
	}
    }
    Base = fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, Base, fsm_parse_regex("[?:0]^4 [?:0 ?:0 ? ?]* [?:0]^4", NULL, NULL))));
    Base = fsm_unflatten(b_handle, Base, "@0@", "@ID@");

    for (i = 0; specialsymbols[i] != NULL; i++) {
	Base->sigma = sigma_remove(specialsymbols[i], Base->sigma);
    }
    for (rule_number = 1; rule_number <= num_rules; rule_number++)
	Base->sigma = sigma_remove(rb->namestrings[rule_number-1], Base->sigma);

    fsm_compact(b_handle, Base);
    sigma_sort(Base);
    rewrite_cleanup(rb);
    return Base;
}

void rewrite_cleanup(struct rewrite_batch *rb) {

    if (rb->Rulenames != NULL)
	fsm_destroy(rb->Rulenames);
    if (rb->ISyms != NULL)
	fsm_destroy(rb->ISyms);
    if (rb->ANY != NULL)
	fsm_destroy(rb->ANY);
    if (rb->IOpen != NULL)
	fsm_destroy(rb->IOpen);
    if (rb->IClose != NULL)
	fsm_destroy(rb->IClose);
    if (rb->ITape != NULL)
	fsm_destroy(rb->ITape);
    if (rb->Any4Tape != NULL)
	fsm_destroy(rb->Any4Tape);
    if (rb->Epextend != NULL)
	fsm_destroy(rb->Epextend);
    if (rb->namestrings != NULL)
	xxfree(rb->namestrings);
    xxfree(rb);
    return;
}


struct fsm *rewr_notlongest(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number, int arrow_type) {
    /* define NotLongest(X)  [Upper(X)/Lower(X) & Tape1of4(IOpen Tape1Sig* ["@O@" | IOpen] Tape1Sig*)] */
    struct fsm *nl, *flt, *rulenum;
    nl = fsm_parse_regex("[\"@I[@\"|\"@I[]@\"] [\"@I[@\"|\"@I[]@\"|\"@I]@\"|\"@I@\"|\"@O@\"]* [\"@O@\"|\"@I[@\"|\"@I[]@\"] [\"@I[@\"|\"@I[]@\"|\"@I]@\"|\"@I@\"|\"@O@\"]*", NULL, NULL);
    nl = rewrite_tape_m_to_n_of_k(b_handle, nl, 1, 1, 4);
    rulenum = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_symbol(rb->namestrings[rule_number-1]), fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_identity(), fsm_universal())))));
    nl = fsm_intersect(b_handle, nl, rulenum);
    /* lang can't end in @0@ */
    if (arrow_type & ARROW_RIGHT) {
	flt = fsm_parse_regex("[? ? ? ?]* [? ? [?-\"@0@\"] ?]", NULL, NULL);
    } else {
	flt = fsm_parse_regex("[? ? ? ?]* [? ? ? [?-\"@0@\"]]", NULL, NULL);
    }
    return fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_intersect(b_handle, nl, fsm_copy(lang)), flt));
}

struct fsm *rewr_notshortest(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number) {
    /* define NotShortest(X)   [Upper/Lower(X) & Tape1of4("@I[@" \IClose*)] */
    struct fsm *ns, *rulenum;
    ns = fsm_parse_regex("[\"@I[@\"] \\[\"@I]@\"]*", NULL, NULL);
    rulenum = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_symbol(rb->namestrings[rule_number-1]), fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_identity(), fsm_universal())))));
    ns = rewrite_tape_m_to_n_of_k(b_handle, ns, 1, 1, 4);
    ns = fsm_intersect(b_handle, ns, rulenum);
    return fsm_minimize(b_handle, fsm_intersect(b_handle, ns, fsm_copy(lang)));
}

struct fsm *rewr_notleftmost(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rule_number, int arrow_type) {
    struct fsm *nl, *flt, *rulenum;
    /* define Leftmost(X)   [Upper/Lower(X) & Tape1of4("@O@" Tape1Sig* IOpen Tape1Sig*) ] */
    nl = fsm_parse_regex("\"@O@\" [\"@O@\"]* [\"@I[@\"|\"@I[]@\"] [\"@I[@\"|\"@I[]@\"|\"@I]@\"|\"@I@\"|\"@O@\"]*", NULL, NULL);
    nl = rewrite_tape_m_to_n_of_k(b_handle, nl, 1, 1, 4);
    rulenum = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_concat(b_handle, fsm_symbol("@O@"), fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_identity(), fsm_identity()))), fsm_concat(b_handle, fsm_kleene_star(b_handle, fsm_concat(b_handle, fsm_symbol("@O@"), fsm_concat(b_handle, fsm_identity(), fsm_concat(b_handle, fsm_identity(), fsm_identity())))), fsm_concat(b_handle, fsm_union(fsm_symbol("@I[@"), fsm_symbol("@I[]@")), fsm_concat(b_handle, fsm_symbol(rb->namestrings[rule_number-1]), fsm_universal())))));
    nl = fsm_intersect(b_handle, nl, rulenum);
    if (arrow_type & ARROW_RIGHT) {
	flt = fsm_parse_regex("[? ? ? ?]* [? ? [?-\"@0@\"] ?]", NULL, NULL); 
    } else {
	flt = fsm_parse_regex("[? ? ? ?]* [? ? ? [?-\"@0@\"]]", NULL, NULL);
    }
    return fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_intersect(b_handle, nl, fsm_copy(lang)), flt));
}

struct fsm *rewr_unrewritten(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang) {
    /* define Unrewritten(X) [X .o. [0:"@O@" 0:"@0@" ? 0:"@ID@"]*].l; */
    struct fsm *C;
    C = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@O@")), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@0@")), fsm_concat(b_handle, fsm_copy(rb->ANY), fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@ID@")))))));
    return fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, lang, C)));
}

struct fsm *rewr_contains(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang) {
    /* define NotContain(X) ~[[Tape1Sig Tape2Sig Tape3Sig Tape4Sig]* X ?*]; */
    return fsm_minimize(b_handle, fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), fsm_concat(b_handle, lang, rewrite_any_4tape(b_handle, rb))));
}

struct fsm *rewrite_tape_m_to_n_of_k(struct build_handle *b_handle, struct fsm *lang, int m, int n, int k) {
    /* [X .o. [0:?^(m-1) ?^(n-m+1) 0:?^(k-n)]*].l */
    return fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, lang, fsm_kleene_star(b_handle, fsm_concat(b_handle, fsm_concat_n(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_identity()), m-1), fsm_concat(b_handle, fsm_concat_n(b_handle, fsm_identity(), n-m+1), fsm_concat_n(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_identity()), k-n)))))));
}

struct fsm *rewrite_two_level(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lang, int rightside) {
    struct fsm *Lower, *Upper, *Result;
    Lower = rewrite_lower(b_handle, rb, fsm_minimize(b_handle, fsm_lower(b_handle, fsm_copy(lang))));
    Upper = rewrite_upper(b_handle, rb, fsm_minimize(b_handle, fsm_upper(b_handle, lang)));
    if (rightside == 1) {
	Result = fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_concat(b_handle, Lower, rewrite_any_4tape(b_handle, rb)), fsm_concat(b_handle, Upper, rewrite_any_4tape(b_handle, rb))));
    } else {
	Result = fsm_minimize(b_handle, fsm_intersect(b_handle, fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), Lower), fsm_concat(b_handle, rewrite_any_4tape(b_handle, rb), Upper)));
    }
    return Result;
}

struct fsm *rewrite_lower(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *lower) {

    /*
       Lower:

       [ @O@      | ISyms    | ISyms    ]*
       [ @0@      | Rulenums | Rulenums ]
       [ <R>,@#@  | @0@,R    |  R       ]
       [ @ID@     | <R>      | @0@      ]

       R = any real symbol
       <R> = any real symbol, not inserted

    */

    struct fsm *One, *Two, *Three, *Filter;

    One = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@O@")), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@0@")), fsm_concat(b_handle, fsm_union(fsm_symbol("@#@"), fsm_copy(rb->ANY)), fsm_cross_product(b_handle, fsm_empty_string(),fsm_symbol("@ID@"))))));

    Two = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ISyms)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->Rulenames)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_union(fsm_copy(rb->ANY), fsm_symbol("@0@"))), fsm_copy(rb->ANY)))));

    Three = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ISyms)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->Rulenames)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ANY)), fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@0@"))))));

    Filter = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_union(One, fsm_union(Two, Three))));
    return fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, lower, Filter)));
}

struct fsm *rewrite_any_4tape(struct build_handle *b_handle, struct rewrite_batch *rb) {

    /*
      Upper:

      [ @O@      | ISyms      ]*
      [ @0@      | Rulenums   ]
      [ <R>,@#@  | @0@,R      ]
      [ @ID@     | R,@ID@,@0@ ]

      R = any real symbol
      <R> = any real symbol, not inserted
    */
    if (rb->Any4Tape == NULL) {
	rb->Any4Tape = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_union(fsm_concat(b_handle, fsm_symbol("@O@"), fsm_concat(b_handle, fsm_symbol("@0@"), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@#@")), fsm_symbol("@ID@")))), fsm_concat(b_handle, fsm_copy(rb->ISyms), fsm_concat(b_handle, fsm_copy(rb->Rulenames), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@0@")), fsm_union(fsm_copy(rb->ANY), fsm_union(fsm_symbol("@ID@"), fsm_symbol("@0@")))))))));
    }
    return fsm_copy(rb->Any4Tape);
}

struct fsm *rewrite_upper(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper) {
    /*
      Upper:

      [ @O@      | ISyms    | ISyms      ]*
      [ @0@      | Rulenums | Rulenums   ]
      [ <R>,@#@  | @0@      | <R>        ]
      [ @ID@     |  R       | R,@ID@,@0@ ]

      R = any real symbol
      <R> = any real symbol, not inserted
    */

    struct fsm *One, *Two, *Three, *Filter;

    One = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@O@")), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@0@")), fsm_concat(b_handle, fsm_union(fsm_symbol("@#@"), fsm_copy(rb->ANY)), fsm_cross_product(b_handle, fsm_empty_string(),fsm_symbol("@ID@"))))));

    Two = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ISyms)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->Rulenames)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_symbol("@0@")), fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ANY))))));

    Three = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->ISyms)), fsm_concat(b_handle, fsm_cross_product(b_handle, fsm_empty_string(), fsm_copy(rb->Rulenames)), fsm_concat(b_handle, fsm_copy(rb->ANY), fsm_cross_product(b_handle, fsm_empty_string(), fsm_union(fsm_union(fsm_symbol("@0@"), fsm_copy(rb->ANY)), fsm_symbol("@ID@")))))));

    Filter = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_union(One, fsm_union(Two, Three))));
    return fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, upper, Filter)));
}

struct fsm *rewrite_align(struct build_handle *b_handle, struct fsm *upper, struct fsm *lower) {
    struct fsm *first, *second, *third, *align, *align2;
    /* `[[`[[Tape1of2(upper "@0@"*) & Tape2of2(lower "@0@"*) & ~[[? ?]* "@0@" "@0@" [? ?]*]], %@%_IDENTITY%_SYMBOL%_%@,%@UNK%@] .o. [? ?|"@UNK@" "@UNK@":"@ID@"]*].l, %@UNK%@,%@%_IDENTITY%_SYMBOL%_%@] */
    first = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, fsm_concat(b_handle, upper, fsm_kleene_star(b_handle, fsm_symbol("@0@"))), 1, 1, 2));
    second = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, fsm_concat(b_handle, lower, fsm_kleene_star(b_handle, fsm_symbol("@0@"))), 2, 2, 2));
    third = fsm_minimize(b_handle, fsm_parse_regex("~[[? ?]* \"@0@\" \"@0@\" [? ?]*]", NULL, NULL));

    align = fsm_minimize(b_handle, fsm_intersect(b_handle, third, fsm_intersect(b_handle, first, second)));
    align = fsm_minimize(b_handle, fsm_substitute_symbol(b_handle, align, "@_IDENTITY_SYMBOL_@", "@UNK@"));
    align2 = fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, align, fsm_parse_regex("[? ? | \"@UNK@\" \"@UNK@\":\"@ID@\" ]*", NULL, NULL))));
    align2 = fsm_minimize(b_handle, fsm_substitute_symbol(b_handle, align2, "@UNK@", "@_IDENTITY_SYMBOL_@"));
    return align2;
}

struct fsm *rewrite_align_markup(struct build_handle *b_handle, struct fsm *upper, struct fsm *lower1, struct fsm *lower2) {
    struct fsm *first, *second, *third, *fourth, *fifth, *sixth, *align, *align2;
    /* [Tape1of2("@0@"*) & Tape2of2(lower1)] [Tape1of2(upper) & Tape2of2("@ID@"*)] [ Tape1of2(lower1) & Tape2of2("@0@"*)] */
    /* + make sure IDENTITY and UNKNOWN are taken care of */
    first = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, fsm_kleene_star(b_handle, fsm_symbol("@0@")), 1, 1, 2));
    second = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, lower1, 2, 2, 2));
    third = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, upper, 1, 1, 2));
    fourth = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, fsm_kleene_star(b_handle, fsm_symbol("@ID@")), 2, 2, 2));
    fifth = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, fsm_kleene_star(b_handle, fsm_symbol("@0@")), 1, 1, 2));
    sixth = fsm_minimize(b_handle, rewrite_tape_m_to_n_of_k(b_handle, lower2, 2, 2, 2));
    align = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_intersect(b_handle, first, second), fsm_concat(b_handle, fsm_intersect(b_handle, third, fourth), fsm_intersect(b_handle, fifth, sixth))));
    align = fsm_minimize(b_handle, fsm_substitute_symbol(b_handle, align, "@_IDENTITY_SYMBOL_@", "@UNK@"));
    align2 = fsm_minimize(b_handle, fsm_lower(b_handle, fsm_compose(b_handle, align, fsm_parse_regex("[? ? | \"@UNK@\" \"@UNK@\":\"@ID@\" ]*", NULL, NULL))));
    align2 = fsm_minimize(b_handle, fsm_substitute_symbol(b_handle, align2, "@UNK@", "@_IDENTITY_SYMBOL_@"));
    return align2;
}

struct fsm *rewrite_itape(struct rewrite_batch *rb) {
    if (rb->ITape == NULL) {
	rb->ITape = fsm_parse_regex("[\"@I[]@\" ? ? ? | \"@I[@\" ? ? ? [\"@I@\" ? ? ?]* \"@I]@\" ? [?-\"@0@\"] ? ] [\"@I]@\" ? \"@0@\" ?]* | 0"  , NULL, NULL);
    }
    return fsm_copy(rb->ITape);
}

struct fsm *rewrite_cp_markup(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper, struct fsm *lower1, struct fsm *lower2, int rule_number) {
    /* Same as rewrite_cp, could be consolidated */
    struct fsm *Aligned, *threetape, *rulenumtape;
    /* define CP(X,Y) Tape23of3(Align2(X,Y)) & [ "@I[@"  ? ? ["@I@" ? ?]* "@I]@" ? ? | "@I[]@" ? ? | 0 ] */
    Aligned = rewrite_align_markup(b_handle, upper, lower1, lower2);
    Aligned = rewrite_tape_m_to_n_of_k(b_handle, Aligned, 3, 4, 4);
    threetape = fsm_minimize(b_handle, fsm_intersect(b_handle, Aligned, rewrite_itape(rb)));
    rulenumtape = rewrite_tape_m_to_n_of_k(b_handle, fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_symbol(rb->namestrings[rule_number-1]))), 2, 2, 4);
    return fsm_minimize(b_handle, fsm_intersect(b_handle, threetape, rulenumtape));
}

struct fsm *rewrite_cp_transducer(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *t, int rule_number) {
    struct fsm *Aligned, *threetape, *rulenumtape;
    Aligned = fsm_flatten(b_handle, t, fsm_symbol("@0@"));
    Aligned = rewrite_tape_m_to_n_of_k(b_handle, Aligned, 3, 4, 4);
    threetape = fsm_minimize(b_handle, fsm_intersect(b_handle, Aligned, rewrite_itape(rb)));
    rulenumtape = rewrite_tape_m_to_n_of_k(b_handle, fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_symbol(rb->namestrings[rule_number-1]))), 2, 2, 4);
    return fsm_minimize(b_handle, fsm_intersect(b_handle, threetape, rulenumtape));
}

struct fsm *rewrite_cp(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *upper, struct fsm *lower, int rule_number) {
    struct fsm *Aligned, *threetape, *rulenumtape;
    /* define CP(X,Y) Tape23of3(Align2(X,Y)) & [ "@I[@"  ? ? ["@I@" ? ?]* "@I]@" ? ? | "@I[]@" ? ? | 0 ] */
    Aligned = rewrite_align(b_handle, upper, lower);
    Aligned = rewrite_tape_m_to_n_of_k(b_handle, Aligned, 3, 4, 4);
    threetape = fsm_minimize(b_handle, fsm_intersect(b_handle, Aligned, rewrite_itape(rb)));
    rulenumtape = rewrite_tape_m_to_n_of_k(b_handle, fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_symbol(rb->namestrings[rule_number-1]))), 2, 2, 4);
    return fsm_minimize(b_handle, fsm_intersect(b_handle, threetape, rulenumtape));
}

void rewrite_add_special_syms(struct rewrite_batch *rb, struct fsm *net) {
    int i;
    if (net == NULL)
        return;
    sigma_substitute(".#.", "@#@", net->sigma); /* We convert boundaries to our interal rep.                          */
                                                /* This is because sigma merging (fsm_merge_sigma()) is handled       */
                                                /* in a special way for .#., which we don't want here.                */

    for (i = 0; specialsymbols[i] != NULL; i++) {
	if (sigma_find(specialsymbols[i], net->sigma) == -1)
	    sigma_add(specialsymbols[i], net->sigma);
    }
    for (i = 1; i <= rb->num_rules; i++) {
	sigma_add(rb->namestrings[i-1], net->sigma);
    }
    sigma_sort(net);
}


void fsm_clear_contexts(struct fsmcontexts *contexts) {
    struct fsmcontexts *c, *cp;
    for (c = contexts; c != NULL; c = cp) {
	fsm_destroy(c->left);
	fsm_destroy(c->right);
	fsm_destroy(c->cpleft);
	fsm_destroy(c->cpright);
	cp = c->next;
	xxfree(c);
    }
}


struct fsm *rewr_context_restrict(struct build_handle *b_handle, struct rewrite_batch *rb, struct fsm *X, struct fsmcontexts *LR) {

    struct fsm *Var, *Notvar, *UnionL, *UnionP, *Result, *NewX, *Left, *Right;
    struct fsmcontexts *pairs;

    Var = fsm_symbol("@VARX@");
    //Notvar = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_term_negation(fsm_symbol("@VARX@"))));
    Notvar = fsm_minus(b_handle, rewrite_any_4tape(b_handle, rb), fsm_contains(b_handle, fsm_symbol("@VARX@")));
    /* We add the variable symbol to all alphabets to avoid ? matching it */
    /* which would cause extra nondeterminism */

    NewX = fsm_copy(X);
    if (sigma_find("@VARX@", NewX->sigma) == -1) {
	sigma_add("@VARX@", NewX->sigma);
	sigma_sort(NewX);
    }
    UnionP = fsm_empty_set();

    for (pairs = LR; pairs != NULL ; pairs = pairs->next) {
	if (pairs->left == NULL) {
	    Left = fsm_empty_string();
	} else {
	    Left = fsm_copy(pairs->cpleft);
	    sigma_add("@VARX@", Left->sigma);
	    sigma_sort(Left);
	}
	if (pairs->right == NULL) {
	    Right = fsm_empty_string();
	} else {
	    Right = fsm_copy(pairs->cpright);
	    sigma_add("@VARX@", Right->sigma);
	    sigma_sort(Right);
	}
        UnionP = fsm_union(fsm_concat(b_handle, Left, fsm_concat(b_handle, fsm_copy(Var), fsm_concat(b_handle, fsm_copy(Notvar), fsm_concat(b_handle, fsm_copy(Var), Right)))), UnionP);
    }
    UnionL = fsm_concat(b_handle, fsm_copy(Notvar), fsm_concat(b_handle, fsm_copy(Var), fsm_concat(b_handle, fsm_copy(NewX), fsm_concat(b_handle, fsm_copy(Var), fsm_copy(Notvar)))));
    Result = fsm_minus(b_handle, UnionL, fsm_concat(b_handle, fsm_copy(Notvar), fsm_concat(b_handle, fsm_copy(UnionP), fsm_copy(Notvar))));

    if (sigma_find("@VARX@", Result->sigma) != -1) {
        Result = fsm_complement(b_handle, fsm_substitute_symbol(b_handle, Result, "@VARX@","@_EPSILON_SYMBOL_@"));
    } else {
        Result = fsm_complement(b_handle, Result);
    }
    fsm_destroy(UnionP);
    fsm_destroy(Var);
    fsm_destroy(Notvar);
    fsm_destroy(NewX);
    return(Result);
}

struct fsm *rewrite_epextend(struct build_handle *b_handle, struct rewrite_batch *rb) {

    struct fsm *one, *two, *allzeroupper, *threea, *threeb, *threec, *three;

    /* 1.  @O@   @0@     [ANY|@#@] @ID@           */
    /* 2.  @I[]@ @#Rule@ [ANY]     [@ID@|@0@|ANY] */
    /* 3a. @I[@  @#Rule@ [ANY]     [@ID@|@0@|ANY] */
    /* 3b. @I@   @#Rule@ [ANY]     [@ID@|@0@|ANY] */
    /* 3c. @I]@  @#Rule@ [ANY]     [@ID@|@0@|ANY] */
    /* 3.  [3a|3b|3c] & ~[[? ? "@0@" ?]*]         */

    /* TODO lower version as well */

    if (rb->Epextend == NULL) {
	one = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_symbol("@O@"), fsm_concat(b_handle, fsm_symbol("@0@"), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@#@")), fsm_symbol("@ID@")))));
	two = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_symbol("@I[]@"), fsm_concat(b_handle, fsm_copy(rb->Rulenames), fsm_concat(b_handle, fsm_copy(rb->ANY), fsm_union(fsm_symbol("@0@"), fsm_union(fsm_symbol("@ID@"), fsm_copy(rb->ANY)))))));
	allzeroupper = fsm_parse_regex("~[[? ? \"@0@\" ?]*]", NULL, NULL);
	threea = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_symbol("@I[@"), fsm_concat(b_handle, fsm_copy(rb->Rulenames), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@0@")), fsm_union(fsm_symbol("@0@"), fsm_union(fsm_symbol("@ID@"), fsm_copy(rb->ANY)))))));
	threeb = fsm_minimize(b_handle, fsm_kleene_star(b_handle, fsm_concat(b_handle, fsm_symbol("@I@"), fsm_concat(b_handle, fsm_copy(rb->Rulenames), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@0@")), fsm_union(fsm_symbol("@0@"), fsm_union(fsm_symbol("@ID@"), fsm_copy(rb->ANY))))))));
	threec = fsm_minimize(b_handle, fsm_concat(b_handle, fsm_symbol("@I]@"), fsm_concat(b_handle, fsm_copy(rb->Rulenames), fsm_concat(b_handle, fsm_union(fsm_copy(rb->ANY), fsm_symbol("@0@")), fsm_union(fsm_symbol("@0@"), fsm_union(fsm_symbol("@ID@"), fsm_copy(rb->ANY)))))));
	three = fsm_intersect(b_handle, allzeroupper, fsm_concat(b_handle, threea, fsm_concat(b_handle, threeb, threec)));
	rb->Epextend = fsm_minimize(b_handle, fsm_union(fsm_union(one, two), three));
    }
    return fsm_copy(rb->Epextend);
}

