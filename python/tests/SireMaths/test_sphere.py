
from Sire.Maths import *

from nose.tools import assert_almost_equal

def _pvt_anaVolume(s0,s1,s2):
    """Analytically calculate the volume of three spheres
       using the inclusion/exclusion principle"""

    v0 = s0.volume()
    v1 = s1.volume()
    v2 = s2.volume()

    v01 = s0.intersectionVolume(s1)
    v02 = s0.intersectionVolume(s2)
    v12 = s1.intersectionVolume(s2)

    v012 = s0.intersectionVolume(s1,s2)

    return v0 + v1 + v2 - v01 - v02 - v12 + v012


def _pvt_checkEqual(val0, val1, verbose):
    delta = abs(val0 - val1)

    if delta > 0.25:
        if verbose:
            print("LARGE DELTA == %s  (%s vs. %s)" % (delta, val0, val1))

        return 1
    else:
        return 0


def test_pair_volume(verbose = False):

    s0 = Sphere(Vector(0), 0.8)
    s1 = Sphere(Vector(-2,0,0), 1.1)

    nerrors = 0

    for i in range(0,41):
        v_an = Sphere.combinedVolume( [s0,s1] )
        v_mc = Sphere.combinedVolumeMC( [s0,s1] )

        if verbose:
            print("%s: %s vs. %s" % (s1.center().x(), v_an, v_mc))

        nerrors += _pvt_checkEqual(v_an, v_mc, verbose)

        s1 = s1.translate( Vector(0.1,0,0) )

    if verbose:
        print("Number of analytic / MC disagreements: %s" % nerrors)

    assert( nerrors < 10 )

def test_triple_volume(verbose = False):

    s0 = Sphere(Vector(0), 0.5)
    s1 = Sphere(Vector(0), 0.7)
    s2 = Sphere(Vector(0), 0.3)

    v_an = Sphere.combinedVolume( [s0,s1,s2] )
    v_ie = _pvt_anaVolume(s0,s1,s2)
    v_mc = Sphere.combinedVolumeMC( [s0,s1,s2], 5000000 )

    if verbose:
        print("%s vs. %s vs. %s" % (v_an, v_ie, v_mc))

    assert_almost_equal( v_an, v_ie, 5 )

    nerrors = 0
    nerrors += _pvt_checkEqual(v_an, v_mc, verbose)

    for i in range(0,25):
        s1 = s1.translate( Vector(0.1,0,0) )

        v_an = Sphere.combinedVolume( [s0,s1,s2] )
        v_ie = _pvt_anaVolume(s0,s1,s2)
        v_mc = Sphere.combinedVolumeMC( [s0,s1,s2], 5000000 )

        if verbose:
            print("%s  %s  %s" % (s0,s1,s2))
            print("%s vs. %s vs. %s" % (v_an, v_ie, v_mc))

        assert_almost_equal(v_an, v_ie, 5)
        nerrors += _pvt_checkEqual(v_an, v_mc, verbose)

    s1 = Sphere( Vector(0.6), 0.7 )

    for i in range(0,25):
        s2 = s2.translate( Vector(0.1,0,0) )

        v_an = Sphere.combinedVolume( [s0,s1,s2] )
        v_ie = _pvt_anaVolume(s0,s1,s2) 
        v_mc = Sphere.combinedVolumeMC( [s0,s1,s2], 5000000 )
    
        if verbose:
            print("%s  %s  %s" % (s0,s1,s2))
            print("%s vs. %s vs. %s" % (v_an, v_ie, v_mc))

        assert_almost_equal(v_an, v_ie, 5)
        nerrors += _pvt_checkEqual(v_an, v_mc, verbose)

    if verbose:
        print("\nNumber of analytic / MC disagreements == %s" % nerrors)

    assert( nerrors < 10 )

if __name__ == "__main__":
    test_pair_volume(True)
    test_triple_volume(True)
