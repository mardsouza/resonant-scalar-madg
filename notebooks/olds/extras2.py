
# pyright: ignore[reportOperatorIssue]

from sympy.combinatorics.permutations import Permutation
import sympy as sy
from ..heppy import *



def vbuwMB(p1, h1, p2, h2, qv, mmed, cv, ca):

    # boson_mass = qv[1]
    qvec = qv / mmed

    # GCW theory
    cv, ca, gz, gw, thetaw, gz_theta, gw_theta, gf = sy.symbols(
        r'c_v c_a g_Z g_W theta_W g_Z_theta g_W_theta, g_f', real=True, )
    # gz = gw / sy.cos(thetaw)
    gv = cv / 2
    ga = ca / 2
    # gz_theta = (gw / sy.cos(thetaw))

    # gw_theta = gw * sy.cos(thetaw)
    EWproj = -sy.I * gf * (gv * one - ga * g5)

    j0 = -sy.I * (vbar(p1, h1) * g0 * EWproj * u(p2, h2) -
               vbar(p1, h1) * g0 * EWproj * qvec[0] * u(p2, h2))
    j1 = -sy.I * (vbar(p1, h1) * g1 * EWproj * u(p2, h2) -
               vbar(p1, h1) * g1 * EWproj * qvec[1] * u(p2, h2))
    j2 = -sy.I * (vbar(p1, h1) * g2 * EWproj * u(p2, h2) -
               vbar(p1, h1) * g2 * EWproj * qvec[2] * u(p2, h2))
    j3 = -sy.I * (vbar(p1, h1) * g3 * EWproj * u(p2, h2) -
               vbar(p1, h1) * g3 * EWproj * qvec[3] * u(p2, h2))
    result = sy.Matrix([sy.simplify(j0[0]), sy.simplify(j1[0]),
                    sy.simplify(j2[0]), sy.simplify(j3[0])])
    return result


def V3ZOutOut(k, hk, p, hp, qv, mmed, cv, ca):

    cv, ca, gz, gw, thetaw, gz_theta, gw_theta, gv = sy.symbols(
        r'c_v c_a g_Z g_W theta_W g_Z_theta g_W_theta, g_v', real=True)

    eps3 = polbar(k, hk)
    eps4 = polbar(p, hp)

    p3 = fourvec(k)
    p4 = fourvec(p)

    # boson_mass = qv[1]
    qvec = qv / mmed

    s1 = dotprod4(eps3, eps4)  # 0
    s2 = dotprod4(eps3, p4)   # GeV
    s3 = dotprod4(eps4, p3)   # GeV

    elm1 = s1 * (p4[0] - p3[0]) - 2 * s2 * eps4[0] + 2 * s3 * eps3[0]  # GeV
    elm2 = s1 * (p4[1] - p3[1]) - 2 * s2 * eps4[1] + 2 * s3 * eps3[1]  # GeV
    elm3 = s1 * (p4[2] - p3[2]) - 2 * s2 * eps4[2] + 2 * s3 * eps3[2]  # GeV
    elm4 = s1 * (p4[3] - p3[3]) - 2 * s2 * eps4[3] + 2 * s3 * eps3[3]  # GeV

    w3 = dotprod4(qvec, p3)     # GeV
    w4 = dotprod4(qvec, p4)     # GeV
    we3 = dotprod4(qvec, eps3)  # 0
    we4 = dotprod4(qvec, eps4)  # 0
    print(f'w3 {w3}')
    print(f'w4 {w4}')
    wd = w3 - w4  # = dotprod4(qvec, p4 - p3)

    SSelm1 = ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
              ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm2 = ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
              ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm3 = ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
              ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm4 = ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
              ((2 * s3) * we4))  # * qvec[0]     ## GeV

    result = sy.Matrix([elm1 - SSelm1, elm2 - SSelm2, elm3 -
                    SSelm3, elm4 - SSelm4]) * sy.I * gv

    return result


def vbuwMB_noqq(p1, h1, p2, h2, qv, mmed, cv, ca):

    # boson_mass = qv[1]
    qvec = qv / mmed

    # GCW theory
    cv, ca, gz, gw, thetaw, gz_theta, gw_theta, gf = sy.symbols(
        r'c_v c_a g_Z g_W theta_W g_Z_theta g_W_theta, g_f', real=True, )
    # gz = gw / sy.cos(thetaw)
    gv = cv / 2
    ga = ca / 2
    # gz_theta = (gw / sy.cos(thetaw))
    # gw_theta = gw * sy.cos(thetaw)
    EWproj = -sy.I * gf * (gv * one - ga * g5)

    # - vbar(p1, h1) * g0 * EWproj * qvec[0]  * u(p2, h2) )
    j0 = -sy.I * (vbar(p1, h1) * g0 * EWproj * u(p2, h2))
    # - vbar(p1, h1) * g1 * EWproj * qvec[1]  * u(p2, h2) )
    j1 = -sy.I * (vbar(p1, h1) * g1 * EWproj * u(p2, h2))
    # - vbar(p1, h1) * g2 * EWproj * qvec[2]  * u(p2, h2) )
    j2 = -sy.I * (vbar(p1, h1) * g2 * EWproj * u(p2, h2))
    # - vbar(p1, h1) * g3 * EWproj * qvec[3]  * u(p2, h2) )
    j3 = -sy.I * (vbar(p1, h1) * g3 * EWproj * u(p2, h2))
    result = sy.Matrix([sy.simplify(j0[0]), sy.simplify(j1[0]),
                    sy.simplify(j2[0]), sy.simplify(j3[0])])
    return result


def V3ZOutOut_noqq(k, hk, p, hp, qv, mmed, cv, ca):

    cv, ca, gz, gw, thetaw, gz_theta, gw_theta, gv = symbols(
        r'c_v c_a g_Z g_W theta_W g_Z_theta g_W_theta, g_v', real=True)

    eps3 = polbar(k, hk)
    eps4 = polbar(p, hp)

    p3 = fourvec(k)
    p4 = fourvec(p)

    qvec = qv / mmed

    s1 = dotprod4(eps3, eps4)  # 0
    s2 = dotprod4(eps3, p4)   # GeV
    s3 = dotprod4(eps4, p3)   # GeV

    elm1 = s1 * (p4[0] - p3[0]) - 2 * s2 * eps4[0] + 2 * s3 * eps3[0]  # GeV  # type: ignore
    elm2 = s1 * (p4[1] - p3[1]) - 2 * s2 * eps4[1] + 2 * s3 * eps3[1]  # GeV  # type: ignore
    elm3 = s1 * (p4[2] - p3[2]) - 2 * s2 * eps4[2] + 2 * s3 * eps3[2]  # GeV  # type: ignore
    elm4 = s1 * (p4[3] - p3[3]) - 2 * s2 * eps4[3] + 2 * s3 * eps3[3]  # GeV  # type: ignore

    w3 = dotprod4(qvec, p3)     # GeV
    w4 = dotprod4(qvec, p4)     # GeV
    we3 = dotprod4(qvec, eps3)  # 0
    we4 = dotprod4(qvec, eps4)  # 0

    wd = w3 - w4  # = dotprod4(qvec, p4 - p3)

    # SSelm1 = ((s1 * (w3 - w4))  - ((2 * s2) * we3)  +  ((2 * s3)  * we4) ) #* qvec[0]     ## GeV
    # SSelm2 = ((s1 * (w3 - w4))  - ((2 * s2) * we3)  +  ((2 * s3)  * we4) ) #* qvec[0]     ## GeV
    # SSelm3 = ((s1 * (w3 - w4))  - ((2 * s2) * we3)  +  ((2 * s3)  * we4) ) #* qvec[0]     ## GeV
    # SSelm4 = ((s1 * (w3 - w4))  - ((2 * s2) * we3)  +  ((2 * s3)  * we4) ) #* qvec[0]     ## GeV

    result = sy.Matrix([elm1, elm2, elm3, elm4]) * sy.I * gv
    # result = sy.Matrix([elm1 - SSelm1, elm2 - SSelm2, elm3 - SSelm3, elm4 - SSelm4]) * sy.I * gv

    return result


# !!!
# Dark matter scenarios

# F + F -> V
def vbu_VDM(p1, h1, p2, h2, qv, mmed, cv, ca, secterm=True):

    st = 1 if secterm else 0

    # boson_mass = qv[1]
    qvec = qv / mmed

    gf, gl, gr, grv, gra, glv, gla = sy.symbols(
        r'g_f g_l g_r g_{rv} g_{ra} g_{lv} g_{la}', real=True, )

    # GCW theory

    # gz = gw / sy.cos(thetaw)
    # glv = cv / 2
    # gla = ca / 2
    # gz_theta = (gw / sy.cos(thetaw))

    # gw_theta = gw * sy.cos(thetaw)
    EWproj = -sy.I * gf * (glv * one - gla * g5)
    Glr = - sy.I * (gl * (glv * one - gla * g5) + gr * (grv * one + gra * g5)) / 2

    # EXPLANATION:
    # gr = 0 -> usual matter, no right handed stuff
    # gl = 0 -> exotic matter, only right handed stuff
    # gl = gr -> mixing, two things coexist
    # gl = - gr -> ???
    # glv, gla, gra, grv vector or axial vector stuff, given a right or left handed scenario, in a right handded scenario replaced by ca and cv, should give EW theory

    j0 = -sy.I * (vbar(p1, h1) * g0 * Glr * u(p2, h2) - st *
               vbar(p1, h1) * g0 * Glr * qvec[0] * u(p2, h2))
    j1 = -sy.I * (vbar(p1, h1) * g1 * Glr * u(p2, h2) - st *
               vbar(p1, h1) * g1 * Glr * qvec[1] * u(p2, h2))
    j2 = -sy.I * (vbar(p1, h1) * g2 * Glr * u(p2, h2) - st *
               vbar(p1, h1) * g2 * Glr * qvec[2] * u(p2, h2))
    j3 = -sy.I * (vbar(p1, h1) * g3 * Glr * u(p2, h2) - st *
               vbar(p1, h1) * g3 * Glr * qvec[3] * u(p2, h2))
    result = sy.Matrix([sy.simplify(j0[0]), sy.simplify(j1[0]),
                    sy.simplify(j2[0]), sy.simplify(j3[0])])
    return result


def V3ZOutOut_VDM(k, hk, p, hp, qv, mmed, cv, ca, secterm=True):

    cv, ca, gz, gw, thetaw, gz_theta, gw_theta, gv = sy.symbols(
        r'c_v c_a g_Z g_W theta_W g_Z_theta g_W_theta, g_v', real=True)

    st = 1 if secterm else 0

    eps3 = polbar(k, hk)
    eps4 = polbar(p, hp)

    p3 = fourvec(k)
    p4 = fourvec(p)

    # boson_mass = qv[1]
    qvec = qv / mmed

    s1 = dotprod4(eps3, eps4)  # 0
    s2 = dotprod4(eps3, p4)   # GeV
    s3 = dotprod4(eps4, p3)   # GeV

    elm1 = s1 * (p4[0] - p3[0]) - 2 * s2 * eps4[0] + 2 * s3 * eps3[0]  # GeV
    elm2 = s1 * (p4[1] - p3[1]) - 2 * s2 * eps4[1] + 2 * s3 * eps3[1]  # GeV
    elm3 = s1 * (p4[2] - p3[2]) - 2 * s2 * eps4[2] + 2 * s3 * eps3[2]  # GeV
    elm4 = s1 * (p4[3] - p3[3]) - 2 * s2 * eps4[3] + 2 * s3 * eps3[3]  # GeV

    w3 = dotprod4(qvec, p3)     # GeV
    w4 = dotprod4(qvec, p4)     # GeV
    we3 = dotprod4(qvec, eps3)  # 0
    we4 = dotprod4(qvec, eps4)  # 0
    print(f'w3 {w3}')
    print(f'w4 {w4}')
    wd = w3 - w4  # = dotprod4(qvec, p4 - p3)

    SSelm1 = st * ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
                   ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm2 = st * ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
                   ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm3 = st * ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
                   ((2 * s3) * we4))  # * qvec[0]     ## GeV
    SSelm4 = st * ((s1 * (w3 - w4)) - ((2 * s2) * we3) +
                   ((2 * s3) * we4))  # * qvec[0]     ## GeV

    result = sy.Matrix([elm1 - SSelm1, elm2 - SSelm2, elm3 -
                    SSelm3, elm4 - SSelm4]) * sy.I * gv

    return result


# F + F -> S
def ff_SDM_(p1, h1, p2, h2, qv, mmed):

    # st = 1 if secterm else 0

    # mediator mass = qv[1]
    qvec = qv / mmed

    gf, gs = sy.symbols(r'g_f g_s', real=True)

    # GCW theory

    # gz = gw / sy.cos(thetaw)
    # glv = cv / 2
    # gla = ca / 2
    # gz_theta = (gw / sy.cos(thetaw))

    # gw_theta = gw * sy.cos(thetaw)
    # EWproj = -sy.I * gf * (glv * one - gla * g5)
    # * (glv * one - gla * g5)  + grx * (grv * one + gra * g5)) / 2
    GS = - sy.I * (gs)

    # EXPLANATION:
    # gr = 0 -> usual matter, no right handed stuff
    # gl = 0 -> exotic matter, only right handed stuff
    # gl = gr -> mixing, two things coexist
    # gl = - gr -> ???
    # glv, gla, gra, grv vector or axial vector stuff, given a right or left handed scenario, in a right handded scenario replaced by ca and cv, should give EW theory

    # * sy.I #- st * vbar(p1, h1) * g0 * Glr * qvec[0]  * u(p2, h2) )
    j0 = (vbar(p1, h1) * GS * u(p2, h2))
    # * sy.I #- st * vbar(p1, h1) * g1 * Glr * qvec[1]  * u(p2, h2) )
    j1 = (vbar(p1, h1) * GS * u(p2, h2))
    # * sy.I #- st * vbar(p1, h1) * g2 * Glr * qvec[2]  * u(p2, h2) )
    j2 = (vbar(p1, h1) * GS * u(p2, h2))
    # * sy.I #- st * vbar(p1, h1) * g3 * Glr * qvec[3]  * u(p2, h2) )
    j3 = (vbar(p1, h1) * GS * u(p2, h2))
    result = sy.Matrix([sy.simplify(j0[0]), sy.simplify(j1[0]),
                    sy.simplify(j2[0]), sy.simplify(j3[0])])
    return result


# S -> F + -F
def _SDM_ff(p1, h1, p2, h2, qv, mmed):

    # st = 1 if secterm else 0

    # mediator mass = qv[1]
    qvec = qv / mmed

    gf, gsx = sy.symbols(r'g_f g_sx', real=True)

    # * (glv * one - gla * g5)  + grx * (grv * one + gra * g5)) / 2
    GS = - sy.I * (gsx)

    # EXPLANATION:
    # gr = 0 -> usual matter, no right handed stuff
    # gl = 0 -> exotic matter, only right handed stuff
    # gl = gr -> mixing, two things coexist
    # gl = - gr -> ???
    # glv, gla, gra, grv vector or axial vector stuff, given a right or left handed scenario, in a right handded scenario replaced by ca and cv, should give EW theory

    j0 = (vbar(p1, h1) * GS * u(p2, h2))
    j1 = (vbar(p1, h1) * GS * u(p2, h2))
    j2 = (vbar(p1, h1) * GS * u(p2, h2))
    j3 = (vbar(p1, h1) * GS * u(p2, h2))
    result = sy.Matrix([sy.simplify(j0[0]), sy.simplify(j1[0]),
                    sy.simplify(j2[0]), sy.simplify(j3[0])])
    return result


def phi(pp, ha):
    res = sy.Matrix([[1], [1], [1], [1]])
    return res


def phibar(pp, ha):
    res = sy.Matrix([[1, 1, 1, 1]])
    return res


def Amplitude_schannel(u1, p1, u2, p2, prp, u3, p3, u4, p4):

    # Breit-Wigner denominator, to be added in the end of the calculation for simplicity.
    spinor = {
        'u': {'type': 'F', 'fn': u},
        'v': {'type': 'F', 'fn': v},
        'vbar': {'type': 'F', 'fn': vbar},
        'ubar': {'type': 'F', 'fn': ubar},
        'pol': {'type': 'V', 'fn': pol},
        'polbar': {'type': 'V', 'fn': polbar},
        'phi': {'type': 'S', 'fn': phi},
        'phibar': {'type': 'S', 'fn': phibar}
    }

    # prop = {
    #         'scalar'    : 'S',
    #         'vector'    : 'V',
    #         'axial'     : 'A',
    #         'pseudoscalar' : 'P',
    #         'fermion' : 'F',
    # }

    # Breit-Wigner denominator, to be added in the end of the calculation for simplicity.
    def get_propagator(prop, type1, dm=False):

        # pick the correct indice for SM or DM mechanics
        dmsm = r"\chi" if dm else r"\psi"

        s, Mmed, Gamma = sy.symbols(r's M_{med} Gamma', real=True)
        gs, gps = sy.symbols(f'g_s{dmsm} g_ps{dmsm}', real=True)
        gv, gav = sy.symbols(f'g_v{dmsm} g_av{dmsm}', real=True)
        # gs, gps, gsx =  symbols(r'g_{v} g_{ps} g_{sx}', real=True )

        denom = (s - Mmed**2 + sy.I * (Mmed * Gamma))

        if prop == 'scalar':
            vertex = {
                # scalar and pseudoscalar couplings
                f"F": -sy.I * (gs * one + sy.I * gps * g5),
                f"V": -sy.I * gs * one,  # Scalar vector coupling
                f"S": -sy.I * gs * one,  # Triple scalar simple coupling
            }

            propagator = -1 / denom

        elif prop == 'vector':
            # vertex = {
            #             f"{type1}{type2}": -sy.I * ( gs * one + sy.I * gps * g5 ), ## scalar and pseudoscalar couplings
            #             f"{type1}{type2}": -sy.I *  gsx * one,
            #             f"{type1}{type2}": -sy.I *  gsx * one,
            #          }                  ## Triple scalar simple coupling
            pass

        elif prop == 'fermion':
            print('Not yet')
            pass

        return vertex[type1], propagator

    def init_J():
        hel_fermion = [-1, 1]
        vertex, propagator = get_propagator(
            prp, spinor[u1]['type'])  # Get initial propagator
        Jc_init = []  # define initial current
        for hl in hel_fermion:
            for hr in hel_fermion:
                jc = []
                for i in range(4):
                    # print(spinor[u1](p1, hl))
                    # print('newspinor')
                    # IMPORTANT!!! intial current - propagator only here
                    ji = spinor[u1]['fn'](
                        p1, hl) * vertex * propagator * spinor[u2]['fn'](p2, hr)
                    jc.append(sy.simplify(ji[0]))

                Jc_init.append(sy.Matrix(jc))

        return Jc_init

    def final_J(types):
        hel_fermion = [-1, 1] if types != 'V' else [-1, 0, 1]
        vertex, propagator = get_propagator(
            prp, types, dm=True)  # Get initial propagator
        Jc_final = []  # define initial current
        for hl in hel_fermion:
            for hr in hel_fermion:
                jc = []
                for i in range(4):
                    # print(spinor[u3](p3, hl))
                    # IMPORTANT!!! # second current
                    jf = spinor[u3]['fn'](p3, hl) * \
                        vertex * spinor[u4]['fn'](p4, hr)
                    # print(jf)
                    jc.append(sy.simplify(jf[0]))

                Jc_final.append(sy.Matrix(jc))

        return Jc_final

    # Trace final result
    Ji = init_J()
    Jf = final_J(spinor[u4]['type'])
    Terms = []
    for j_init in Ji:
        for j_final in Jf:
            # print('append term')

            Terms.append(sy.simplify(dotprod4(j_init, j_final)))

    return Terms


def decay_Gamma(prp, u3, p3, u4, p4):

    # Breit-Wigner denominator, to be added in the end of the calculation for simplicity.
    spinor = {
        'u': {'type': 'F', 'fn': u},
        'v': {'type': 'F', 'fn': v},
        'vbar': {'type': 'F', 'fn': vbar},
        'ubar': {'type': 'F', 'fn': ubar},
        'pol': {'type': 'V', 'fn': pol},
        'polbar': {'type': 'V', 'fn': polbar},
        'phi': {'type': 'S', 'fn': phi},
        'phibar': {'type': 'S', 'fn': phibar}
    }

    # prop = {
    #         'scalar'    : 'S',
    #         'vector'    : 'V',
    #         'axial'     : 'A',
    #         'pseudoscalar' : 'P',
    #         'fermion' : 'F',
    # }

    # Breit-Wigner denominator, to be added in the end of the calculation for simplicity.
    def get_propagator(prop, types, dm=False):

        # pick the correct indice for SM or DM mechanics
        dmsm = r"\chi" if dm else r"\psi"

        s, Mmed, Gamma = sy.symbols(r's M_{med} Gamma', real=True)
        gs, gps = sy.symbols(f'g_s{dmsm} g_ps{dmsm}', real=True)
        gv, gav = sy.symbols(f'g_v{dmsm} g_av{dmsm}', real=True)
        # gs, gps, gsx =  symbols(r'g_{v} g_{ps} g_{sx}', real=True )

        denom = (s - Mmed**2 + sy.I * (Mmed * Gamma))

        if prop == 'scalar':
            vertex = {
                # + sy.I * gps * g5 ), ## scalar and pseudoscalar couplings
                f"F": sy.I * -sy.I * (gs * one),
                f"V": sy.I * gs * one,  # Scalar vector coupling
                f"S": -sy.I * gs * one,  # Triple scalar simple coupling
            }

            propagator = -1 / denom

        elif prop == 'vector':
            # vertex = {
            #             f"{type1}{type2}": -sy.I * ( gs * one + sy.I * gps * g5 ), ## scalar and pseudoscalar couplings
            #             f"{type1}{type2}": -sy.I *  gsx * one,
            #             f"{type1}{type2}": -sy.I *  gsx * one,
            #          }                  ## Triple scalar simple coupling
            pass

        elif prop == 'fermion':
            print('Not yet')
            pass

        return vertex[types], propagator

    def init_J(types):
        # print(types)
        hel_fermion = [-1, 1] if types != 'V' else [-1, 0, 1]
        vertex, propagator = get_propagator(
            prp, types)  # Get initial propagator
        Jc_init = []  # define initial current
        for hl in hel_fermion:
            jc = []
            for i in range(4):
                # print(hl)
                # print(spinor[u1](p1, hl))
                # print('newspinor')
                ji = spinor[u3]['fn'](p3, hl) * vertex
                jc.append(sy.simplify(ji[0]))

            Jc_init.append(sy.Matrix(jc))

        return Jc_init

    def final_J(types):
        # print(types)
        hel_fermion = [-1, 1] if types != 'V' else [-1, 0, 1]
        vertex, propagator = get_propagator(
            prp, types, dm=True)  # Get initial propagator
        Jc_final = []  # define initial current

        for hl in hel_fermion:
            jc = []

            for i in range(4):
                # print(hl)
                # print(spinor[u3](p3, hl))
                jf = one * spinor[u4]['fn'](p4, hl)
                # print(jf)
                jc.append(sy.simplify(jf[0]))

            Jc_final.append(sy.Matrix(jc))

        return Jc_final

    # Trace final result
    Ji = init_J(spinor[u3]['type'])
    Jf = final_J(spinor[u4]['type'])
    Terms = []
    for j_init in Ji:
        for j_final in Jf:
            # print('append term')

            Terms.append(sy.simplify(dotprod4(j_init, j_final)))

    return Terms


