from columnflow.util import maybe_import
from columnflow.columnar_util import attach_coffea_behavior
ak = maybe_import("awkward")

# creates record with useful variables
def table_combo(array: ak.Array) -> ak.Array:
    table = array.metric_table(array, 
                            metric=lambda a, b: abs((a+b).mass-125))
    # create dictionary with indicies
    idx1 = ak.local_index(table, axis=-2)
    idx2 = ak.local_index(table, axis=-1)
    table_combo = ak.zip({"mass_diff": table,
                          "idx1": idx1, "idx2": idx2})
    # remove duplicates and self sums
    table_combo = ak.mask(table_combo,table_combo.idx1<table_combo.idx2)
    return table_combo

# gets all unique index permutations for picking two set of pairs    
def pair_permutations(array: ak.Array) -> ak.Array:
    # # create all unique pair permutaions
    idx = ak.local_index(array, axis=1)
    pairs = ak.combinations(idx, 2, fields=["idx1","idx2"],axis=1)
    permu = ak.combinations(pairs, 2,fields=["pair1","pair2"],axis=1)
    permu_mask = ((permu.pair1.idx1 != permu.pair2.idx1) & 
                  (permu.pair1.idx1 != permu.pair2.idx2) & 
                  (permu.pair1.idx2 != permu.pair2.idx1) & 
                  (permu.pair1.idx2 != permu.pair2.idx2))
    permu = permu[permu_mask]
    permu = ak.pad_none(permu,1)
    return permu

# puts out the two pairs where chi**2 is minimized and also the value of chi**2
def min_chi_sqr_pair(array: ak.Array, table: ak.Array) -> ak.Array:
    # minimise chi**2 for pairings
    permu = pair_permutations(array)
    chisq = ((((array[permu.pair1.idx1] + array[permu.pair1.idx2]).mass -125)/125)**2 + 
             (((array[permu.pair2.idx1] + array[permu.pair2.idx2]).mass - 125)/125)**2)
    sorted_chi_idx = ak.argsort(chisq, axis=1, ascending=True)
    bestpairs = permu[sorted_chi_idx][:,0]
    pairs_mask = (((table.idx1 == bestpairs.pair1.idx1) & (table.idx2 == bestpairs.pair1.idx2)) | 
                  ((table.idx1 == bestpairs.pair2.idx1) & (table.idx2 == bestpairs.pair2.idx2)))
    chi_table = ak.mask(table, pairs_mask)
    chi_table = ak.flatten(ak.drop_none(chi_table,axis=1),axis=2)

    return chi_table
