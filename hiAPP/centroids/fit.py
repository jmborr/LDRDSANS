from mantid import simpleapi as msapi
import scipy.cluster.hierarchy as sch
import navigate

def fit(sim_i, sim_q, exp_i, exp_q):
    """
    Fit simulated I(Q) against experimental
    :param sim_i: (numpy.array) simulated profile
    :param sim_q: (numpy.array) Q-values where sim_i is evaluated
    :param exp_i: (numpy.array) experimental profile
    :param exp_q: (numpy.array) Q-values where exp_i is evaluated
    :return: chi-square value
    """
    simw = msapi.CreateWorkspace(DataX=sim_q, DataY=sim_i, UnitX='MomentumTransfer')
    expw = msapi.CreateWorkspace(DataX=exp_q, DataY=exp_i, UnitX='MomentumTransfer')
    function = 'TabulatedFunction(Workspace="simw")+FlatBackground'
    return Fit(Function=function, InputWorkspace='expw')

if __name__ == '__main__':
    """
    Fit average I(Q) of each cluster against experimental I(Q)
    """
    centroids_I = msapi.LoadNexus("centroids.sans.nxs")  # Load I(Q) for each centroid
    Z = np.loadtxt('./centroids.linkage_matrix')
    root_node = sch.to_tree(Z)
    for cs_at_depth in navigate.clusters_afloat_depth(root_node, depth=10):
        for c in cs_at_depth:
        sim_i, sim_q = calculate_average_I(c.)
            results = fit(sim_i, sim_q, exp_i, exp_q)
            fits.append(((depth, id_at_depth), results))