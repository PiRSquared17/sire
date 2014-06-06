
from Sire.Maths import *
from Sire.Units import *

from math import sin,cos

from nose.tools import assert_almost_equal

def _pvt_makeRotateX(angle):
    rad = angle.to(radians)

    return Matrix( 1, 0, 0,
                   0, cos(rad), -sin(rad),
                   0, sin(rad), cos(rad) )

def _pvt_makeRotateY(angle):
    rad = angle.to(radians)

    return Matrix(cos(rad), 0, sin(rad),
                  0, 1, 0,
                  -sin(rad), 0, cos(rad) )

def _pvt_makeRotateZ(angle):
    rad = angle.to(radians)

    return Matrix(cos(rad), -sin(rad), 0,
                  sin(rad), cos(rad), 0,
                  0, 0, 1)

def _assert_vector_equal(v0, v1, n=3):
    assert_almost_equal( v0.x(), v1.x(), n )
    assert_almost_equal( v0.y(), v1.y(), n )
    assert_almost_equal( v0.z(), v1.z(), n )

def _assert_matrix_equal(m0, m1, n=3):
    for i in range(0,3):
        for j in range(0,3):
            assert_almost_equal( m0.at(i,j), m1.at(i,j), n )

def _assert_quaternion_equal(q0, q1, n=3):

    try:
        assert_almost_equal( q0.x(), q1.x(), n )
        assert_almost_equal( q0.y(), q1.y(), n )
        assert_almost_equal( q0.z(), q1.z(), n )
        assert_almost_equal( q0.w(), q1.w(), n )
    except:
        try:
            assert_almost_equal( q0.x(), -q1.x(), n )
            assert_almost_equal( q0.y(), -q1.y(), n )
            assert_almost_equal( q0.z(), -q1.z(), n )
            assert_almost_equal( q0.w(), -q1.w(), n )
        except:
            print("Are these quaternions equal? %s vs. %s" % (q0,q1))
            _assert_matrix_equal(q0.toMatrix(), q1.toMatrix())

x = Vector(1,0,0)
y = Vector(0,1,0)
z = Vector(0,0,1)
w = Vector(1,1,1)

def _test_rotate(func, v, verbose):
    for i in range(0, 361):
        m = func( i * degrees )
        q = Quaternion( i * degrees, v )

        if verbose:
            print("Testing rotation %s degrees about %s..." % (i,v))
            print("%s vs. %s" % (Quaternion(m), q))
            print(m)
            print("vs.")
            print(q.toMatrix())

        _assert_quaternion_equal( q, Quaternion(q.toMatrix()) )
        _assert_matrix_equal( m, Quaternion(m).toMatrix() )
        _assert_matrix_equal( m, q.toMatrix() )
        _assert_quaternion_equal( q, Quaternion(m) )

        if verbose:
            print("Rotate x: %s vs. %s" % (m*x, q.rotate(x)))
            print("Rotate y: %s vs. %s" % (m*y, q.rotate(y)))
            print("Rotate z: %s vs. %s" % (m*z, q.rotate(z)))
            print("Rotate w: %s vs. %s" % (m*w, q.rotate(w)))

        _assert_vector_equal( m * x, q.rotate(x) )
        _assert_vector_equal( m * y, q.rotate(y) )
        _assert_vector_equal( m * z, q.rotate(z) )    
        _assert_vector_equal( m * w, q.rotate(w) )    

def _test_rotate_q(v, verbose):
    for i in range(0, 361):
        q = Quaternion( i * degrees, v )
        m = q.toMatrix()

        if verbose:
            print("Testing rotation %s degrees about %s..." % (i,v))
            print("%s vs. %s" % (Quaternion(m), q))
            print(m)
            print("vs.")
            print(q.toMatrix())

        if i != 180:
            # pathological error
            _assert_quaternion_equal( q, Quaternion(q.toMatrix()) )
            _assert_matrix_equal( m, Quaternion(m).toMatrix() )
            _assert_matrix_equal( m, q.toMatrix() )
            _assert_quaternion_equal( q, Quaternion(m) )

        if verbose:
            print("Rotate x: %s vs. %s" % (m*x, q.rotate(x)))
            print("Rotate y: %s vs. %s" % (m*y, q.rotate(y)))
            print("Rotate z: %s vs. %s" % (m*z, q.rotate(z)))
            print("Rotate w: %s vs. %s" % (m*w, q.rotate(w)))

        _assert_vector_equal( m * x, q.rotate(x) )
        _assert_vector_equal( m * y, q.rotate(y) )
        _assert_vector_equal( m * z, q.rotate(z) )        
        _assert_vector_equal( m * w, q.rotate(w) )    

def test_rotate(verbose = False):

    # test all rotations around x, y, and z
    _test_rotate(_pvt_makeRotateX, x, verbose)
    _test_rotate(_pvt_makeRotateY, y, verbose)
    _test_rotate(_pvt_makeRotateZ, z, verbose)

    # now try some more complex rotations
    _test_rotate_q( Vector(1,1,1), verbose )
    

if __name__ == "__main__":
    test_rotate(True)
