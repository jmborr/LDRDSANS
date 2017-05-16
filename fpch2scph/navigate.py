import numpy as np

def clusters_afloat_depth(node, depth=-1):
    """
    Return list of clusters at or less than depth for node
    :param node: ClusterNode
    :param depth: maximum depth from node. Default value 0f -1 indicates depth = node.count-1
    :return: list of lists, each sublist contains ClusterNode objects at depths less
     or equal than depth. The first list is always [node,]
    """
    maximum_depth = node.count - 1  # all clusters are leafs at maximum_depth
    if depth < 0 or depth > maximum_depth:
        depth = maximum_depth
    clusters_at_all_depths = [[node,]]
    for d in range(1, depth+1):
        clusters_at_d = list(clusters_at_all_depths[d-1])
        max_id = np.argmax([cluster.id for cluster in clusters_at_d])
        c = clusters_at_d[max_id]  # split the largest cluster
        clusters_at_d.extend((c.left, c.right))
        clusters_at_d.remove(c)
        clusters_at_all_depths.append(clusters_at_d)
    return clusters_at_all_depths

def leafs_under_node(node):
    """
    List of leafs under node
    :param node: ClusterNode
    :return: list of ClusterNode leafs
    """
    return node.pre_order(lambda x: x)


if __name__ == '__main__':
    import scipy.cluster.hierarchy as sch
    import pytest
    np.random.seed(4711)  # for repeatability of this tutorial
    a = np.random.multivariate_normal([10, 0], [[3, 1], [1, 4]], size=[100, ])
    b = np.random.multivariate_normal([0, 20], [[3, 1], [1, 4]], size=[50, ])
    X = np.concatenate((a, b), )
    Z = sch.linkage(X, 'ward')
    root_node = sch.to_tree(Z)
    cad = clusters_afloat_depth(root_node)
    leafs = leafs_under_node(root_node)

