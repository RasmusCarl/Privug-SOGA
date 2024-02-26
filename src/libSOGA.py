# Contains the function start_soga, which is used to invoke SOGA on a CFG object and the recursive function SOGA, which, depending on the type of the visited node, calls the functions needed to update the current distribution. 

#Such functions are contained in the auxiliary libraries:
# - libSOGAtruncate, containing functions for computing the resulting distribution when a truncation occurs (in conditional or observe instructions);
# - libSOGAupdate, containing functions for computing the resulting distribution after applying an assignment instruction;
# - libSOGAmerge, containing functions for computing the resulting distribution when a merge instruction is encountered;

# Additional functions for general purpose are defined in the library libSOGAshared, which is imported by all previous libraries.

# TO DO:
# - improve dependencies on libraries (all auxiliary libraries import libSOGAshared, maybe there is a more efficient way to do this?)
# - libSOGAmerge: add other pruning  

from libSOGAtruncate import *
from libSOGAupdate import *
from libSOGAmerge import *
import timing

def start_SOGA(cfg, pruning=None, Kmax=None, parallel=False,useR=False):
    """ Invokes SOGA on the root of the CFG object cfg, initializing current_distribution to a Dirac delta centered in zero.
        If pruning='classic' implements pruning at the merge nodes with maximum number of component Kmax.
        Returns an object Dist (defined in libSOGAshared) with the final computed distribution."""
    #timing.trunc_time.clear()
    #timing.update_time.clear()
    #timing.prune_time.clear()
    #timing.num_trunc = 0
    if(useR):
        initR()

    var_list = cfg.ID_list
    data = cfg.data
    gm = GaussianMix([1.], [np.array([0.]*len(var_list))], [np.zeros((len(var_list),len(var_list)))])
    init_dist = Dist(var_list, gm)
    SOGA(cfg.root, 1, init_dist, None, data, parallel)
    p, current_dist = merge(cfg.node_list['exit'].list_dist)
    #print('Truncations:', timing.num_trunc, 'total:', sum(timing.trunc_time))
    #print('Updates:', len(timing.update_time), 'total:', sum(timing.update_time))
    #print('Pruning:', len(timing.prune_time), 'total:', sum(timing.prune_time))
    cfg.node_list['exit'].list_dist = []
    return current_dist


def SOGA(node, current_p, current_dist, current_trunc, data, parallel):
    
    #print('Current #components:', current_dist.gm.n_comp())
        
    # starts execution
    if node.type == 'entry':
        for child in node.children:
            SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
    
    # if tests saves LBC and calls on children
    if node.type == 'test':
        current_trunc = node.LBC
        for child in node.children:
            SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
            
    # if loop saves checks the condition and decides which child node must be accessed
    if node.type == 'loop':
        # the first time is accessed set the value of the counter to 0 and converts node.const into a number
        if data[node.idx][0] is None:
            data[node.idx][0] = 0
        if type(node.const) is str:
            if '[' in node.const:
                data_name, data_idx = node.const.split('[')
                data_idx = data_idx[:-1]
                # data_idx is a data
                if data_idx in data:
                    data_idx = data[data_idx][0]
                # data_idx is a number
                else:
                    data_idx = int(data_idx)
                node.const = int(data[data_name][data_idx])
            else:
                node.const = int(node.const)            
        # successively checks the condition and decides which child node must be accessed
        if data[node.idx][0] < node.const:
            for child in node.children:
                if child.cond == True:
                    SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
        else:
            data[node.idx][0] = None
            for child in node.children:
                if child.cond == False:
                    SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
     
    # if state checks wheter cond!=None. If yes, truncates to current_trunc, eventually negating it. In any case applies the rule in expr. Appends the distribution in the next merge node or calls recursively on children. If child is loop node increments its idx.
    if node.type == 'state':
        if node.cond != None and not current_trunc is None:
            if node.cond == False:
                current_trunc = negate(current_trunc) 
            timing.num_trunc += current_dist.gm.n_comp()
            start = time()
            if parallel:
                p, current_dist = parallel_truncate(current_dist, current_trunc, data)   ### see libSOGAtruncate
            else:
                p, current_dist = truncate(current_dist, current_trunc, data) 
            current_trunc = None
            end = time()
            timing.trunc_time.append(end-start)
            current_p = p*current_p
        if current_p > prob_tol:
            start = time()
            current_dist = update_rule(current_dist, node.expr, data)         ### see libSOGAupdate
            end = time()
            timing.update_time.append(end-start)
        if node.children[0].type == 'merge' or node.children[0].type == 'exit':
            node.children[0].list_dist.append((current_p, current_dist))
        if node.children[0].type == 'loop' and not data[node.children[0].idx][0] is None:
            data[node.children[0].idx][0] += 1
        SOGA(node.children[0], deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
            
    # if observe truncates to LBC and calls on children
    if node.type == 'observe':
        current_trunc = node.LBC
        start = time()
        if parallel:
            p, current_dist = parallel_truncate(current_dist, current_trunc, data)
        else:
            p, current_dist = truncate(current_dist, current_trunc, data)
        #current_p = current_p*p
        current_trunc = None
        end = time()
        timing.trunc_time.append(end-start)
        for child in node.children:
            if child.type == 'merge' or child.type == 'exit':
                child.list_dist.append((current_p, current_dist))
            SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)

    # if merge checks whether all paths have been explored.
    # Either returns or merge distributions and calls on children
    if node.type == 'merge':
        if len(node.list_dist) != len(node.parent):
            return
        else:
            current_p, current_dist = merge(node.list_dist)        ### see libSOGAmerge
            node.list_dist = []
            for child in node.children:
                if child.type == 'merge' or child.type == 'exit':
                    child.list_dist.append((current_p, current_dist))
                SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)
    if node.type == 'exit':
        return
    
    if node.type == 'prune':
        start = time()
        current_dist = prune(current_dist,'classic',node.Kmax)        ### see libSOGAmerge
        end = time()
        timing.prune_time.append(end-start)
        node.list_dist = []
        for child in node.children:
            if child.type == 'merge' or child.type == 'exit':
                child.list_dist.append((current_p, current_dist))
            SOGA(child, deepcopy(current_p), deepcopy(current_dist), deepcopy(current_trunc), data, parallel)

