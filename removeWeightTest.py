import unittest
#trend u_s >> u_a

def test_removeWeight():
    assertequal(removeWeight(1, u_a = 3,  u_s = 400), 400/403)
    assertequal(removeWeight(0.33, u_a = 0.2,  u_s = 76.6), 0.329141)
    #assuming input1 is weight, 2 is absorption coefficient, & 3 is scattering coefficient
    return 