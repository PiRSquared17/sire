
#include "expression.h"
#include "exbase.h"
#include "sum.h"
#include "product.h"
#include "symbol.h"
#include "symbols.h"
#include "functions.h"
#include "power.h"
#include "powerconstant.h"
#include "complexvalues.h"
#include "i.h"
#include "integrationconstant.h"

#include "SireMaths/maths.h"
#include "SireMaths/complex.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;
using namespace SireMaths;

static const RegisterMetaType<Expression> r_expression;

/** Serialise an Expression to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Expression &ex)
{
    writeHeader(ds, r_expression, 1) << ex.fac << ex.exbase;

    return ds;
}

/** Deserialise an Expression from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Expression &ex)
{
    VersionID v = readHeader(ds, r_expression);

    if (v == 1)
    {
        ds >> ex.fac >> ex.exbase;
    }
    else
        throw version_error(v, "1", r_expression, CODELOC);

    return ds;
}

/** Construct an empty (zero) expression */
Expression::Expression() : fac(0)
{}

/** Construct a constant expression equal to 'constant' */
Expression::Expression(int constant) : fac(constant)
{}

/** Construct a constant expression equal to 'constant' */
Expression::Expression(double constant) : fac(constant)
{}

/** Construct a constant expression equal to 'constant' */
Expression::Expression(const Rational &constant) : fac(SireMaths::toDouble(constant))
{}

/** Construct a constant expression equal to 'constant' */
Expression::Expression(const Complex &constant)
{
    if (constant.isZero())
        fac = 0.0;
    else if (constant.isReal())
        fac = constant.real();
    else if (SireMaths::isZero(constant.real()))
    {
        fac = constant.imag();
        exbase = I().toExpressionBase();
    }
    else
    {
        fac = constant.imag();
        exbase = Sum((constant.real() / fac), I());
    }
}

/** Construct an expression that is equal to 1*(base) */
Expression::Expression(const ExpressionBase &base)
           : exbase(base), fac(1)
{
    if (exbase.isConstant() and not exbase.isComplex())
    {
        exbase = ExpressionBase();
        fac = base.evaluate(Values());
    }
}

/** Construct an expression that is equal to 1*(base) */
Expression::Expression(const ExBase &base)
           : exbase( base.toExpressionBase() ), fac(1)
{
    if (exbase.isConstant() and not exbase.isComplex())
    {
        exbase = ExpressionBase();
        fac = base.evaluate(Values());
    }
}

/** Copy constructor */
Expression::Expression(const Expression &other)
           : exbase(other.exbase), fac(other.fac)
{}

/** Destructor */
Expression::~Expression()
{}

/** Assignment operator */
Expression& Expression::operator=(const Expression &other)
{
    exbase = other.exbase;
    fac = other.fac;
    return *this;
}

/** self-addition */
Expression& Expression::operator+=(const Expression &other)
{
    return operator=( this->add(other) );
}

/** self-subtraction */
Expression& Expression::operator-=(const Expression &other)
{
    return operator=( this->subtract(other) );
}

/** self-multiplication */
Expression& Expression::operator*=(const Expression &other)
{
    return operator=( this->multiply(other) );
}

/** self-division */
Expression& Expression::operator/=(const Expression &other)
{
    return operator=( this->divide(other) );
}

/** Evaluate the numerical value of this expression, using the values
    supplied in 'values'. Any unidentified symbols or functions are
    assumed to be equal to zero. Note that this only performs real-arithmetic,
    so an exception will be thrown if any part of this expression generates
    a complex result.

    \throw SireMaths::domain_error
*/
double Expression::evaluate(const Values &values) const
{
    return fac * exbase.evaluate(values);
}

/** Evaluate the numerical value of this expression using complex
    arithmetic. Any unidentified symbols or functions are assumed
    to be equal to zero.
*/
Complex Expression::evaluate(const ComplexValues &values) const
{
    return fac * exbase.evaluate(values);
}

/** Convienient synonym for evaluate() */
double Expression::operator()(const Values &values) const
{
    return evaluate(values);
}

/** Convienient synonym for evaluate() */
Complex Expression::operator()(const ComplexValues &values) const
{
    return evaluate(values);
}

/** Return this expression raised to the power 'n' */
Expression Expression::pow(int n) const
{
    if (n == 1)
        return *this;
    else if (n == 0)
        return Expression(1);
    else if (this->isConstant())
        return SireMaths::pow(evaluate(ComplexValues()), n);
    else
        return IntegerPower( *this, n ).reduce();
}

/** Return this expression raised to the rational power 'n' */
Expression Expression::pow(const Rational &n) const
{
    if (n.denominator() == 1)
        return this->pow(n.numerator());
    else if (this->isConstant())
        return SireMaths::pow(evaluate(ComplexValues()), n);
    else
        return RationalPower( *this, n ).reduce();
}

/** Return this expression raised to a real number power */
Expression Expression::pow(double n) const
{
    if (SireMaths::areEqual(n,1.0))
        return *this;
    else if (this->isConstant())
        return SireMaths::pow(evaluate(ComplexValues()), n);
    else
        return RealPower( *this, n ).reduce();
}

/** Return this expresssion raised to a complex power */
Expression Expression::pow(const Complex &n) const
{
    if (n.isReal())
        return pow(n.real());
    else if (this->isConstant())
        return SireMaths::pow(evaluate(ComplexValues()), n);
    else
        return ComplexPower( *this, n ).reduce();
}

/** Return this expression raised to a function */
Expression Expression::pow(const Expression &n) const
{
    if (n.isConstant())
        return pow(n.evaluate(ComplexValues()));
    else
        return Power( *this, n ).reduce();
}

/** Return 1 / expression */
Expression Expression::invert() const
{
    return this->pow(-1);
}

/** Return the square of this expression */
Expression Expression::squared() const
{
    return this->pow(2);
}

/** Return the cube of this expression */
Expression Expression::cubed() const
{
    return this->pow(3);
}

/** Return the nth root of this expression */
Expression Expression::root(int n) const
{
    return this->pow( Rational(1,n) );
}

/** Return the negative of this expression */
Expression Expression::negate() const
{
    Expression ret(*this);
    ret.fac = -(ret.fac);
    return ret;
}

/** Negation operator */
Expression Expression::operator-() const
{
    return this->negate();
}

/** Return this expression added to 'ex' */
Expression Expression::add(const Expression &ex) const
{
    return Sum(*this, ex).reduce();
}

/** Return an expression that is this - ex */
Expression Expression::subtract(const Expression &ex) const
{
    return Sum(*this, ex.negate()).reduce();
}

/** Return an expression that is this multipled by 'val' */
Expression Expression::multiply(double val) const
{
    if (exbase.isA<IntegrationConstant>())
    {
        return Expression(exbase);
    }
    else
    {
        double newfactor = fac * val;

        if (SireMaths::isZero(newfactor))
            return Expression(0);
        else
        {
            Expression ret(*this);
            ret.fac = newfactor;

            return ret;
        }
    }
}

/** Return an expression that is this multiplied by the complex value z */
Expression Expression::multiply(const Complex &z) const
{
    if (z.isReal())
        return multiply(z.real());
    else if (SireMaths::isZero(z.real()))
        return Product( *this, z.imag() * I() ).reduce();
    else
        return Product( *this, z.real() + z.imag()*I() ).reduce();
}

/** Return an expression that is this multiplied by 'ex' */
Expression Expression::multiply(const Expression &ex) const
{
    if (this->isConstant())
        return ex.multiply( this->evaluate(ComplexValues()) );
    else if (ex.isConstant())
        return multiply( ex.evaluate(ComplexValues()) );
    else
        return Product(*this, ex).reduce();
}

/** Return an expression that is this divided by 'val' */
Expression Expression::divide(double val) const
{
    if ( SireMaths::areEqual(val,1.0) )
        return *this;
    else
        return multiply( double(1.0) / val );
}

/** Return an expression that is divided by the complex number z */
Expression Expression::divide(const Complex &z) const
{
    return multiply( z.inverse() );
}

/** Return an expression that is this / ex */
Expression Expression::divide(const Expression &ex) const
{
    if (ex.isConstant())
        return divide( ex.evaluate(ComplexValues()) );
    else
        return multiply(ex.invert());
}

/** Return an expression whereby the identities in 'identities' have
    been substituted into this expression */
Expression Expression::substitute(const Identities &identities) const
{
    Expression ret = exbase.substitute(identities);
    ret.fac *= fac;
    return ret;
}

/** Return an expanded form of this expression (where possible products of
    sums multiplied out to form sums of products) */
Expression Expression::expand() const
{
    Sum sum(*this);
    return sum.expand();
}

/** Return a collapsed form of this expression (common factors collected together,
    where possible sums of products collected together into products of sums) */
Expression Expression::collapse() const
{
    Product product(*this);
    return product.collapse();
}

/** Return the complex conjugate of this expression */
Expression Expression::conjugate() const
{
    return fac * exbase.conjugate();
}

/** Try to simplify this expression by using built-in identities. If
    SireCAS::UNSAFE_COMPLEX_SIMPLIFICATIONS is passed, then allow the use
    of identities that are not necessarily true in the complex domain,
    e.g. z = sin(arcsin(z)) */
Expression Expression::simplify(int options) const
{
    return fac * exbase.simplify(options);
}

/** Differentiate this expression with respect to 'symbol' and return
    the resulting expression.

    \throw SireCAS::unavailable_differential
*/
Expression Expression::differentiate(const Symbol &symbol, int level) const
{
    if (level <= 0 or this->isZero())
        return *this;
    else if (level > 0 and this->isConstant())
        return Expression(0);
    else if (not exbase.isFunction(symbol))
        // d f(y) / dx = 0
        return Expression(0);
    else
    {
        //calculate the differential of the base expression with respect
        //to symbol
        Expression diff = fac * exbase.differentiate(symbol);

        if (not diff.isZero())
        {
            if (level == 1)
                return diff;
            else
                return diff.differentiate(symbol, level-1);
        }
        else
            return diff;
    }
}

/** Synonym for differentiate */
Expression Expression::diff(const Symbol &symbol, int level) const
{
    return differentiate(symbol,level);
}

/** Integrate this expression with respect to 'symbol' and return the
    resulting expression.

    \throw SireCAS::unavailable_integral
*/
Expression Expression::integrate(const Symbol &symbol) const
{
    if (this->isZero())
        return Expression(0);
    else if (not exbase.isFunction(symbol))
        //exbase is constant with respect to 'symbol' - return symbol*exbase + C
        return fac * exbase + IntegrationConstant();
    else
        //calculate the integral with respect to the symbol (add integration constant)
        return fac * exbase.integrate(symbol) + IntegrationConstant();
}

/** Synonym for integrate */
Expression Expression::integ(const Symbol &symbol) const
{
    return integrate(symbol);
}

/** Return a series expansion of this function with respect to 'symbol'
    up to order 'n'. If an expansion is not possible, then this just
    returns this expression */
Expression Expression::series(const Symbol &symbol, int n) const
{
    return fac * exbase.series(symbol,n);
}

/** Return whether or not this expression is equal to '0' for all values */
bool Expression::isZero() const
{
    return SireMaths::isZero(fac);
}

/** Return whether or not this expression is constant for all values */
bool Expression::isConstant() const
{
    return exbase.isConstant();
}

/** Return whether or not this is a function of 'symbol' */
bool Expression::isFunction(const Symbol &symbol) const
{
    return exbase.isFunction(symbol);
}

/** Return whether or not this is a compound expression (contains more
    than a single expression, e.g. Sum, Product or PowerFunction) */
bool Expression::isCompound() const
{
    return exbase.isCompound();
}

/** Return whether or not this expression has complex parts */
bool Expression::isComplex() const
{
    return exbase.isComplex();
}

/** Return a string representation of this expression */
QString Expression::toString() const
{
    if (isConstant())
        return evaluate(ComplexValues()).toString();
    else
    {
        QString basestr = exbase.toString();

        if ( SireMaths::areEqual(fac,1.0) )
            return basestr;
        else if ( SireMaths::areEqual(fac,-1.0) )
            return QString("-%1").arg(basestr);
        else if (this->isCompound())
            return QString("%1 [%2]").arg(fac).arg(basestr);
        else
            return QString("%1 %2").arg(fac).arg(basestr);
    }
}

/** Return the ExpressionBase base-part of this expression */
const ExpressionBase& Expression::base() const
{
    return exbase;
}

/** Return the factor of this expression */
double Expression::factor() const
{
    return fac;
}

/** Return the child expressions that make up this expression */
Expressions Expression::children() const
{
    return exbase.children();
}

/** Return all of the symbols used in this expression */
Symbols Expression::symbols() const
{
    return exbase.symbols();
}

/** Return all of the functions used in this expression */
Functions Expression::functions() const
{
    return exbase.functions();
}
