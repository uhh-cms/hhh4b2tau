from columnflow.util import maybe_import
from columnflow.columnar_util import attach_coffea_behavior
ak = maybe_import("awkward")
np = maybe_import("numpy")

# creates record with useful variables
def table_combo(array: ak.Array, idx: ak.Array) -> ak.Array:
    array = array[idx]
    table = array.metric_table(
        array, 
        metric=lambda a, b: (((a+b).mass-125)/125)**2
        )
    # create dictionary with indicies
    idx1 = idx
    idx2 = ak.unflatten(idx, 1)[ak.zeros_like(idx)]
    table_combo = ak.zip({"chi2": table,
                          "idx1": idx1, "idx2": idx2})
    # remove duplicates and self sums
    table_combo = ak.mask(table_combo,table_combo.idx1<table_combo.idx2)
    return table_combo

# gets all unique index permutations for picking two set of pairs    
def pair_permutations(array: ak.Array, idx: ak.Array) -> ak.Array:
    # # create all unique pair permutaions
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
def min_chi_sqr_pair(array: ak.Array, table: ak.Array, idx: ak.Array) -> ak.Array:
    # minimise chi**2 for pairings
    permu = pair_permutations(array, idx)
    chisq = chi2(array, permu)
    sorted_chi_idx = ak.argsort(chisq, axis=1, ascending=True)
    best_pairs = permu[sorted_chi_idx][:,0]
    pairs_mask = (
        ((table.idx1 == best_pairs.pair1.idx1) & (table.idx2 == best_pairs.pair1.idx2)) | 
        ((table.idx1 == best_pairs.pair2.idx1) & (table.idx2 == best_pairs.pair2.idx2)) |
        ((table.idx1 == best_pairs.pair1.idx2) & (table.idx2 == best_pairs.pair1.idx1)) | 
        ((table.idx1 == best_pairs.pair2.idx2) & (table.idx2 == best_pairs.pair2.idx1))
        )
    chi_table = ak.mask(table, pairs_mask)
    chi_table = ak.flatten(ak.drop_none(chi_table,axis=1),axis=2)

    min_chisq = chisq[sorted_chi_idx][:,0]
    # rando_mask = (np.random.rand(len(array)) >= 0.5)
    # bb1 = ak.where(rando_mask, chi_table[:,0], chi_table[:,1])
    # bb2 = ak.where(rando_mask, chi_table[:,1], chi_table[:,0])

    bb1 = chi_table[:,0]
    bb2 = chi_table[:,1]

    return bb1, bb2, min_chisq

def mds(array: ak.Array, permu: ak.Array) -> ak.Array:
    mds_1 = (((array[permu.pair1.idx1] + array[permu.pair1.idx2]).mass -125)/125)**2
    mds_2 = (((array[permu.pair2.idx1] + array[permu.pair2.idx2]).mass -125)/125)**2
    mds = ak.concatenate([mds_1, mds_2], axis=1)
    return mds


def chi2(array: ak.Array, permu: ak.Array) -> ak.Array:
    chi2 = ((((array[permu.pair1.idx1] + array[permu.pair1.idx2]).mass -125)/125)**2 + 
             (((array[permu.pair2.idx1] + array[permu.pair2.idx2]).mass - 125)/125)**2)
    return chi2


def dhh(array: ak.Array, permu: ak.Array) -> ak.Array:
    # reconstruction function taken from HHH --> 4b2g analysis
    # ratio of mean masses between H1 and H2
    a = 1.05 
    # D_{HH} is to be minimized
    dhh1 = abs(
        (array[permu.pair1.idx1] + array[permu.pair1.idx2]).mass 
        - a * (array[permu.pair2.idx1] + array[permu.pair2.idx2]).mass
        ) / (1 + a**2)
    dhh2 = abs(
        (array[permu.pair2.idx1] + array[permu.pair2.idx2]).mass 
        - a * (array[permu.pair1.idx1] + array[permu.pair1.idx2]).mass
        ) / (1 + a**2)
    
    dhh = ak.concatenate([dhh1,dhh2],axis=1)
    
    return dhh

def min_func_pair(func, array: ak.Array, table: ak.Array, idx: ak.Array, ordered = False) -> ak.Array:
    # minimise func for pairings
    # ordered == True for functions that have intrinsic ordering so the indicies exceed object numbers
    permu = pair_permutations(array, idx)
    func_val = func(array, permu)
    min_func_val_idx = ak.argsort(func_val, axis=1, ascending=True)[:, :1]

    if ordered == True:
        overload_mask = (min_func_val_idx >= ak.num(permu))
        overload_idx = (min_func_val_idx - ak.num(permu))
        min_func_val_idx = ak.where(overload_mask, overload_idx, min_func_val_idx)


    best_pairs = ak.firsts(permu[min_func_val_idx])

    pair1_mask = (
        ((table.idx1 == best_pairs.pair1.idx1) & (table.idx2 == best_pairs.pair1.idx2)) | 
        ((table.idx1 == best_pairs.pair1.idx2) & (table.idx2 == best_pairs.pair1.idx1))
    )

    pair2_mask = (
        ((table.idx1 == best_pairs.pair2.idx1) & (table.idx2 == best_pairs.pair2.idx2)) | 
        ((table.idx1 == best_pairs.pair2.idx2) & (table.idx2 == best_pairs.pair2.idx1))
    )
    
    pair1 = ak.mask(table, pair1_mask)
    pair1 = ak.flatten(ak.drop_none(pair1, axis=1), axis=2)

    pair2 = ak.mask(table, pair2_mask)
    pair2 = ak.flatten(ak.drop_none(pair2, axis=1), axis=2)

    if ordered == True:
        bb1, bb2 = order_pairs(pair1, pair2, ~overload_mask)

    elif ordered == False:
        bb1, bb2 =swap_random(pair1, pair2)

    min_func_val = ak.min(func_val, axis=1)

    bb1_idx = get_idx(bb1)
    bb2_idx = get_idx(bb2)

    return bb1_idx, bb2_idx, min_func_val


# randomize bb1 and bb2
def swap_random(array0: ak.Array, array1: ak.Array) -> ak.Array:
    rando_mask = (np.random.rand(len(array0)) >= 0.5)
    bb1 = ak.where(rando_mask, array0, array1)
    bb2 = ak.where(rando_mask, array1, array0)
    return bb1, bb2

def get_idx(array: ak.Array) -> ak.Array:
    idx = ak.concatenate([array.idx1, array.idx2],axis=1)
    idx = ak.drop_none(idx)
    return idx

def order_pairs(array0: ak.Array, array1: ak.Array, mask: ak.Array) -> ak.Array:
    pair1 = ak.where(mask, array0, array1)
    pair2 = ak.where(mask, array1, array0)
    return pair1, pair2
