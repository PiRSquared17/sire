/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2009  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/

#include "gto.h"

#include "SireMaths/maths.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

using namespace Squire;
using namespace SireMaths;
using namespace SireBase;
using namespace SireID;
using namespace SireStream;

//////////
////////// Implementation of GTO
//////////

static const RegisterMetaType<GTO> r_gto( MAGIC_ONLY, GTO::typeName() );

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const GTO &gto)
{
    writeHeader(ds, r_gto, 1);
    
    ds << gto.alfa << gto.scl << static_cast<const OrbitalShell&>(gto);
    
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, GTO &gto)
{
    VersionID v = readHeader(ds, r_gto);
    
    if (v == 1)
    {
        ds >> gto.alfa >> gto.scl >> static_cast<OrbitalShell&>(gto);
    }
    else
        throw version_error(v, "1", r_gto, CODELOC);
        
    return ds;
}

/** Constructor */
GTO::GTO() : OrbitalShell(), alfa(0), scl(0)
{}

/** Construct with the specified value of the alpha (exponent)
    and scale (unnormalised) */
GTO::GTO(double alpha, double scale) 
    : OrbitalShell(), alfa(alpha), scl(scale)
{}

/** Copy constructor */
GTO::GTO(const GTO &other) 
    : OrbitalShell(other), alfa(other.alfa), scl(other.scl)
{}

/** Destructor */
GTO::~GTO()
{}

/** Copy assignment operator */
GTO& GTO::operator=(const GTO &other)
{
    if (this != &other)
    {
        alfa = other.alfa;
        scl = other.scl;
        OrbitalShell::operator=(other);
    }
    
    return *this;
}

/** Comparison operator */
bool GTO::operator==(const GTO &other) const
{
    return alfa == other.alfa and scl == other.scl and
           OrbitalShell::operator==(other);
}

/** Comparison operator */
bool GTO::operator!=(const GTO &other) const
{
    return not this->operator==(other);
}

/** Return the value of alpha (the exponent) for this gaussian */
double GTO::alpha() const
{
    return alfa;
}

/** Convenient synonym for 'alpha()' - so you can write
    a.alpha() * b.beta() and have it mean what you expect */
double GTO::beta() const
{
    return GTO::alpha();
}

/** Return the unnormalised scaling factor for this gaussian */
double GTO::scale() const
{
    return scl;
}

const char* GTO::typeName()
{
    return "Squire::GTO";
}

//////////
////////// Implementation of ShellPair
//////////

static const RegisterMetaType<GTOPair> r_gtopair( MAGIC_ONLY, GTOPair::typeName() );

/** Serialise to a binary datastream */
QDataStream SQUIRE_EXPORT &operator<<(QDataStream &ds, const GTOPair &gtopair)
{
    writeHeader(ds, r_gtopair, 1);
    
    ds << gtopair._P << gtopair._R2 << gtopair._zeta << gtopair._xi
       << gtopair._K << gtopair._ss
       << static_cast<const ShellPair&>(gtopair);
       
    return ds;
}

/** Extract from a binary datastream */
QDataStream SQUIRE_EXPORT &operator>>(QDataStream &ds, GTOPair &gtopair)
{
    VersionID v = readHeader(ds, r_gtopair);
    
    if (v == 1)
    {
        ds >> gtopair._P >> gtopair._R2 >> gtopair._zeta
           >> gtopair._xi >> gtopair._K >> gtopair._ss
           >> static_cast<ShellPair&>(gtopair);
    }
    else
        throw version_error(v, "1", r_gtopair, CODELOC);
        
    return ds;
}

/** Null constructor */
GTOPair::GTOPair() 
        : ShellPair(), _R2(0), _zeta(0), _xi(0), _K(0), _ss(0)
{}

static const double sqrt_two_times_pi_to_5_4 = std::sqrt(2.0) * std::pow(pi, (5.0/4.0));
static const double four_over_pi2 = 4 / (pi*pi);

/** Construct a shell pair from the passed two GTO orbital shells,
    located at the specified points */
GTOPair::GTOPair(const Vector &A, const GTO &a,
                 const Vector &B, const GTO &b)
        : ShellPair()
{
    //the product of two Gaussians is a Gaussian :-)
    const double alpha_times_beta = a.alpha() * b.alpha();
    const double alpha_plus_beta = a.alpha() + b.alpha();
    
    _P = (a.alpha()*A + b.alpha()*B) / alpha_plus_beta;
    _R2 = Vector::distance2(A, B);
    _zeta = alpha_plus_beta;
    _xi = alpha_times_beta / alpha_plus_beta;
    
    const double scl = a.scale() * b.scale() *
                          std::pow( four_over_pi2 * a.alpha() * b.beta(), 0.75 );
    
    const double G = scl * std::exp(-_xi*_R2);
    _ss = std::pow( pi / _zeta, 1.5 ) * G;

    _K = sqrt_two_times_pi_to_5_4 * G / alpha_plus_beta;
}

/** Copy constructor */
GTOPair::GTOPair(const GTOPair &other) 
        : ShellPair(other),
          _P(other._P), _R2(other._R2), _zeta(other._zeta), _xi(other._xi), 
          _K(other._K), _ss(other._ss)
{}

/** Destructor */
GTOPair::~GTOPair()
{}

const char* GTOPair::typeName()
{
    return "Squire::GTOPair";
}

/** Copy assignment operator */
GTOPair& GTOPair::operator=(const GTOPair &other)
{
    _P = other._P;
    _R2 = other._R2;
    _zeta = other._zeta;
    _xi = other._xi;
    _K = other._K;
    _ss = other._ss;
    
    ShellPair::operator=(other);
    
    return *this;
}

/** Comparison operator */
bool GTOPair::operator==(const GTOPair &other) const
{
    return _P == other._P and _R2 == other._R2 and 
           _zeta == other._zeta and _K == other._K and
           ShellPair::operator==(other);
}

/** Comparison operator */
bool GTOPair::operator!=(const GTOPair &other) const
{
    return not this->operator==(other);
}

/** Return the rho value for the two passed GTOPair pairs 

    rho = (zeta * eta) / (zeta + eta)
*/
double GTOPair::rho(const GTOPair &P, const GTOPair &Q)
{
    return (P.zeta() * Q.eta()) / (P.zeta() + Q.eta());
}

/** Return the T value for the two passed GTOPair pairs

    T = rho (P-Q)^2
*/
double GTOPair::T(const GTOPair &P, const GTOPair &Q)
{
    return GTOPair::rho(P,Q) * Vector::distance2(P.P(), Q.Q());
}

/** Return the prefactor value for the two passed GTOPair pairs

    preFac = K_AB K_CD / Sqrt(zeta+eta)
*/
double GTOPair::preFac(const GTOPair &P, const GTOPair &Q)
{
    return P.K_AB() * Q.K_CD() / std::sqrt(P.zeta() + Q.eta());
}

/** Return the W value for the two passed GTOPair pairs 

    W = (zeta/(zeta+eta)) P + (eta/(zeta+eta)) Q
*/
Vector GTOPair::W(const GTOPair &P, const GTOPair &Q)
{
    const double zeta_plus_eta = P.zeta() + Q.eta();
    
    return ( (P.zeta()/zeta_plus_eta) * P.P() ) + 
           ( (Q.eta()/zeta_plus_eta) * Q.Q() );
}