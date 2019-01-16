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
#include "foma/foma_r.h"

struct fsm *fsm_lower(struct build_handle *b_handle, struct fsm *net) {
    struct fsm_state *fsm;
    int i, prevstate, out;
    fsm = net->states;
    fsm_state_init(b_handle->da_handle, sigma_max(net->sigma));
    prevstate = -1;
    for (i = 0; (fsm+i)->state_no != - 1; prevstate = (fsm+i)->state_no, i++) {
        if (prevstate != -1 && prevstate != (fsm+i)->state_no) {
            fsm_state_end_state(b_handle->da_handle);
        }
        if (prevstate != (fsm+i)->state_no) {
            fsm_state_set_current_state(b_handle->da_handle, (fsm+i)->state_no, (fsm+i)->final_state, (fsm+i)->start_state);
        }
        if ((fsm+i)->target != -1) {
            out = ((fsm+i)->out == UNKNOWN) ? IDENTITY : (fsm+i)->out;
            fsm_state_add_arc(b_handle->da_handle, (fsm+i)->state_no, out, out, (fsm+i)->target, (fsm+i)->final_state, (fsm+i)->start_state);
        }
    }
    fsm_state_end_state(b_handle->da_handle);
    xxfree(net->states);
    fsm_state_close(b_handle->da_handle, net);
    fsm_update_flags(net,NO,NO,NO,UNK,UNK,UNK);
    sigma_cleanup(net,0);
    return(net);
}

struct fsm *fsm_upper(struct build_handle *b_handle, struct fsm *net) {
    struct fsm_state *fsm;
    int i, prevstate, in;
    fsm = net->states;
    fsm_state_init(b_handle->da_handle, sigma_max(net->sigma));
    prevstate = -1;
    for (i = 0; (fsm+i)->state_no != - 1; prevstate = (fsm+i)->state_no, i++) {
        if (prevstate != -1 && prevstate != (fsm+i)->state_no) {
            fsm_state_end_state(b_handle->da_handle);
        }
        if (prevstate != (fsm+i)->state_no) {
            fsm_state_set_current_state(b_handle->da_handle, (fsm+i)->state_no, (fsm+i)->final_state, (fsm+i)->start_state);
        }
        if ((fsm+i)->target != -1) {
            in = ((fsm+i)->in == UNKNOWN) ? IDENTITY : (fsm+i)->in;
            fsm_state_add_arc(b_handle->da_handle, (fsm+i)->state_no, in, in, (fsm+i)->target, (fsm+i)->final_state, (fsm+i)->start_state);
        }
    }
    fsm_state_end_state(b_handle->da_handle);
    xxfree(net->states);
    fsm_state_close(b_handle->da_handle, net);
    fsm_update_flags(net,NO,NO,NO,UNK,UNK,UNK);
    sigma_cleanup(net,0);
    return(net);
}

