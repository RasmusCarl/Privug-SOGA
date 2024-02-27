# Contains the functions for computing the resulting distribution when a merge is invoked according to the following dependencies.

# SOGA (defined in SOGA.py)
# |- merge
#     |- prune
#        |- classic_prune
#            |- dist
#            |- merge_comp


from src.core.libSOGAshared import *


def merge(list_dist):
    """
    Given a list of couples (p,dist), where each dist is a GaussianMix object, computes a couple (current_p, current_dist), in which current_pi is the sum of p and current_dist is a single GaussianMix object.
    """
    final_pi = []
    final_mu = []
    final_sigma = []
    current_p = 0
    # creates a single mixture
    for (p, dist) in list_dist:
        if p > 0:
            current_p += p
            final_pi = final_pi + list(p*np.array(dist.gm.pi))
            final_mu = final_mu + list(dist.gm.mu)
            final_sigma = final_sigma + list(dist.gm.sigma)
    if len(final_pi) == 0:
        d = len(list_dist[0][1].gm.mu[0])
        #print('no components found')
        return 0, Dist(list_dist[0][1].var_list, GaussianMix([0], [np.array([0]*d)], [np.zeros((d,d))]))
    final_pi = list(np.array(final_pi)/current_p)
    # deletes components with probability less than tol
    zero_list = np.where(np.array(final_pi) < prob_tol)[0]
    if len(zero_list)>0:
        if len(zero_list) == len(final_pi):
            d = len(dist.gm.mu[0])
            #print('no components found')
            return 0, Dist(list_dist[0][1].var_list, GaussianMix([0], [np.array([0]*d)], [np.zeros((d,d))]))
        else:
            for index in sorted(zero_list, reverse=True):
                del final_pi[index]
                del final_mu[index]
                del final_sigma[index]
    current_dist = Dist(list_dist[0][1].var_list, GaussianMix(final_pi, final_mu, final_sigma))
    return current_p, current_dist


def prune(current_dist, pruning, Kmax):
    if pruning == 'classic':
        current_dist = classic_prune(current_dist, Kmax)
    return current_dist


def classic_prune(current_dist, Kmax):
    #print('Pruning ', current_dist.gm.n_comp(), ' components to ', Kmax)
    if current_dist.gm.n_comp() > Kmax:
        n = current_dist.gm.n_comp()
        matrixcost = np.ones((n,n))*np.inf
        for i in range(n):
            for j in range(i):
                new_mu = (current_dist.gm.pi[i]*current_dist.gm.mu[i] + current_dist.gm.pi[j]*current_dist.gm.mu[j]) / (current_dist.gm.pi[i]+current_dist.gm.pi[j])
                matrixcost[i,j] = current_dist.gm.pi[i]*dist(current_dist.gm.mu[i],new_mu) + current_dist.gm.pi[j]*dist(current_dist.gm.mu[j],new_mu)
        while n > Kmax:
            min_idx = np.where(matrixcost == np.min(matrixcost))
            i, j = min_idx
            i = i[0]
            j = j[0]
            min_mu = (current_dist.gm.pi[i]*current_dist.gm.mu[i] + current_dist.gm.pi[j]*current_dist.gm.mu[j]) / (current_dist.gm.pi[i]+current_dist.gm.pi[j])
            current_dist = merge_comp(current_dist, i, j, min_mu)
            n = current_dist.gm.n_comp()
            matrixcost = matrix_delete(matrixcost, i, j)
            if n > Kmax:
                for j in range(n-1):
                    new_mu = (current_dist.gm.pi[n-1]*current_dist.gm.mu[n-1] + current_dist.gm.pi[j]*current_dist.gm.mu[j]) / (current_dist.gm.pi[n-1]+current_dist.gm.pi[j])
                    matrixcost[n-1,j] = current_dist.gm.pi[n-1]*dist(current_dist.gm.mu[n-1],new_mu) + current_dist.gm.pi[j]*dist(current_dist.gm.mu[j],new_mu)
    return current_dist

def matrix_delete(matrix, i, j):
    matrix = np.delete(matrix, max(i,j), axis=0)
    matrix = np.delete(matrix, max(i,j), axis=1)
    matrix = np.delete(matrix, min(i,j), axis=0)
    matrix = np.delete(matrix, min(i,j), axis=1)
    matrix = np.vstack([matrix, np.ones((1, len(matrix)))*np.inf])
    matrix = np.hstack([matrix, np.ones((matrix.shape[0],1))*np.inf])
    return matrix

def dist(vec1, vec2):
    return sum((np.array(vec1)-np.array(vec2))**2)


def merge_comp(current_dist, i, j, mu):
    new_pi = []
    new_mu = []
    new_sigma = []
    for k in range(current_dist.gm.n_comp()):
        if k != i and k != j:
            new_pi.append(current_dist.gm.pi[k])
            new_mu.append(current_dist.gm.mu[k])
            new_sigma.append(current_dist.gm.sigma[k])
    pi = current_dist.gm.pi[i] + current_dist.gm.pi[j]
    new_pi.append(pi)
    new_mu.append(mu)
    # creates a GaussianMix with two component i,j and renormalized weights
    double_pi = np.array([current_dist.gm.pi[i], current_dist.gm.pi[j]])/pi
    new_gm = GaussianMix(double_pi, [current_dist.gm.mu[i], current_dist.gm.mu[j]], [current_dist.gm.sigma[i], current_dist.gm.sigma[j]])
    sigma = new_gm.cov()
    new_sigma.append(sigma) 
    return Dist(current_dist.var_list, GaussianMix(new_pi, new_mu, new_sigma))


