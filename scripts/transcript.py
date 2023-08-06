# n: number of nodes, gLen: length of group element
def yurek_dxr_coeff(n, t, gLen):
    return ((t+1)*gLen + n*64)/2**10 

def yurek_dxr_eval(n, t, gLen):
    return (n*gLen + n*64)/2**10 

# Our transcript size with multisignature scheme
def our_msig(n, t, gLen):
    return (n*gLen + n/8 + BLS_SIG_LEN + 64*t)/2**10

# Our transcript size without multisignature scheme
def our_plain_sig(n, t, gLen):
    return (n*gLen + n/8 + (n-t)*ED_SIG_LEN + 64*t)/2**10

def our_msig_sync(n, gLen, f):
    return (n*gLen + n/8 + BLS_SIG_LEN + 64*f)/2**10

def our_plain_sig_sync(n, gLen, f):
    return (n*gLen + n/8 + (n-f)*ED_SIG_LEN + 64*f)/2**10


def class_group(n, gLen):
    return (219+gLen)*n/2**10

def ve_vss(n, g1Len, g2Len, t):
    NUM_CHUNKS = 16
    NUM_ZK_REPETITIONS = 32

    comm_len = n*g1Len
    ctxt_len = (n+1)*g1Len # ElGamal ciphertext of h^r
    sh_pf_len = 2*g1Len + g2Len + 2*32

    cnk_ctxt_len = (NUM_CHUNKS+1)*n*g1Len
    cnk_pf_len = (2*NUM_ZK_REPETITIONS + n + 3)*g1Len

    return (comm_len + ctxt_len + sh_pf_len + cnk_ctxt_len + cnk_pf_len)/2**10

def mixed_vss(n, g1Len, g2Len, t, deg, enc):
    NUM_CHUNKS = 16
    NUM_ZK_REPETITIONS = 32

    comm_len = n*g1Len
    ctxt_len = (enc+1)*g1Len # ElGamal ciphertext of h^r
    sh_pf_len = 2*g1Len + g2Len + 2*32

    cnk_ctxt_len = (NUM_CHUNKS+1)*enc*g1Len
    cnk_pf_len = (2*NUM_ZK_REPETITIONS + enc + 3)*g1Len

    ve_part = comm_len + ctxt_len + sh_pf_len + cnk_ctxt_len + cnk_pf_len
    reveal_count = (n-(2*t+1+enc))
    reveal_part = 64*reveal_count

    sig_part = 32 + n 

    return (ve_part + reveal_part + sig_part)/2**10

ts = [85, 170, 341]
sigLen = 48
gLen = 48
g2Len = 96
BLS_SIG_LEN = 48
ED_SIG_LEN = 64

for t in ts:
    n = 3*t+1
    deg = 2*t
    # print("Yurek coeff:  ", yurek_dxr_coeff(n, t, gLen))
    # print("Yurek eval:  ", yurek_dxr_eval(n, t, gLen))
    # print("Our Plain:  ", our_plain_sig(n, t, gLen))
    # print("Our Msig:   ", our_msig(n, t, gLen))
    
    # print("Our Plain n/2", our_msig_sync(n, gLen, n/2))
    # print("Our Plain n/4", our_msig_sync(n, gLen, n/4))

    # print("Our Msig n/2", our_plain_sig_sync(n, gLen, n/2))
    # print("Our Msig n/4", our_plain_sig_sync(n, gLen, n/4))
    # print(f"VE VSS for n={n}, size={ve_vss(n, gLen, g2Len, t)} KBytes")

    # print("Mixed VSS with multisig:")
    # print(f"For n={n}, deg={deg}, size={mixed_vss(n, gLen, g2Len, t, deg, t)} KBytes")

    print("Mixed VSS with multisig:")
    print(f"For n={n}, size={class_group(n, gLen)} KBytes")

    print("-"*50)