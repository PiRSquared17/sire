
#include "testunits.h"

#include "SireUnits/units.h"
#include "SireUnits/convert.h"

using namespace SireUnits;
using namespace SireTest;

TestUnits::TestUnits()
{}

TestUnits::~TestUnits()
{}

void TestUnits::initialise(test_suite *test)
{
    test->add( BOOST_TEST_CASE( &TestUnits::runTests ) );
}

void TestUnits::runTests()
{
    qDebug() << QObject::tr("Running TestUnits tests...");

    //print out the values of the physical constants
    qDebug() << QObject::tr("PI = %1, 1/4*Pi*Epsilon0 = %2 (UnitCharge-2 A kcal mol-1)")
                    .arg(pi).arg(one_over_four_pi_eps0);
    qDebug() << QObject::tr("Gas constant = %1 (kcal mol-1 K-1)")
                    .arg(gasr);
    qDebug() << QObject::tr("One atmosphere = %1 (kcal mol-1 A-3)").arg(atm);

    //test conversion of epsilon_0 from SI units to internal units
    BOOST_CHECK( testEqual( convertFrom(8.854187817e-12, (coulomb*coulomb) / (joule*meter)), epsilon0 ) );
    
    //test conversion of gas constant from SI units to internal units
    BOOST_CHECK( testEqual( convertFrom(8.31447215, joule / (mole*kelvin)), gasr ) );

    //test conversion of atmospheres into internal pressure units (kcal mol-1 A-3)
    BOOST_CHECK( testEqual( convertFrom(101325.0, joule / (meter*meter*meter)), atm ) );

    //now check that the constant agree with those calculated and used in ProtoMS2
    BOOST_CHECK( testEqual( one_over_four_pi_eps0, 332.06371 ) );
    BOOST_CHECK( testEqual( gasr, 0.00198720657, 1e-4) );
    BOOST_CHECK( testEqual( atm, 1.45839756e-05 ) );

    //check conversions
    BOOST_CHECK( testEqual( convertFrom(5.5, kcal_per_mol), convertFrom(5500.0, cal_per_mol) ) );
    BOOST_CHECK( testEqual( convertFrom(15.31, angstrom), convertFrom(1.531, nanometer) ) );
    
    
}
