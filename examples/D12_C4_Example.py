from stallings_amalgams.groups import create_D2n, create_Cn
from stallings_amalgams import get_red_precover
from stallings_amalgams.visualization import visualize, datatotex

A = {6: 2}

H = [
    ["b"],
    ["a", "b", "a^-1"],
    ["a", "c", "a", "c", "a", "c"],
]

D12 = create_D2n(6)

C4 = create_Cn(4)


P = get_red_precover(D12, C4, A, H)

visualize(P)
