
#include "qhash_sirecas.h"

#include <QMutex>
#include <QHash>

#include "symbol.h"
#include "symbols.h"
#include "functions.h"
#include "expressions.h"
#include "identities.h"
#include "registerexpression.h"
#include "complexvalues.h"
#include "values.h"

#include "SireStream/datastream.h"

using namespace SireStream;
using namespace SireCAS;

/** Register the 'Symbol' class */
static RegisterExpression<Symbol> RegisterSymbol;

typedef struct
{
QMutex idmutex;
QHash<QString,SymbolID> name2id;
SymbolID lastid;
} SymbolReg;

static SymbolReg *registry = 0;

/** Return an ID for the symbol with representation 'rep'. This 
    creates a new ID if there is no symbol currently registered with 
    this ID. */
SymbolID Symbol::getNewID(const QString &rep)
{
    if (rep.isNull() or rep.isEmpty())
        return 0;

    if (registry == 0)
    {
        registry = new SymbolReg();
        registry->lastid = 0;
    }
    
    QMutexLocker lkr(&(registry->idmutex));
    
    if (registry->name2id.contains(rep))
        return registry->name2id.value(rep);
    else
    {
        registry->lastid++;
        registry->name2id.insert( rep, registry->lastid );
        return registry->lastid;
    }
}

static const RegisterMetaType<Symbol> r_symbol("SireCAS::Symbol");

/** Hash a symbol */
uint Symbol::hash() const
{
    return (r_symbol.magicID() << 16) | (id & 0x0000FFFF);
}

/** Serialise a Symbol to a binary datastream */
QDataStream SIRECAS_EXPORT &operator<<(QDataStream &ds, const Symbol &sym)
{
    writeHeader(ds, r_symbol, 1) << sym.stringrep;
    
    return ds;
}

/** Deserialise a Symbol from a binary datastream */
QDataStream SIRECAS_EXPORT &operator>>(QDataStream &ds, Symbol &sym)
{
    VersionID v = readHeader(ds, r_symbol);
    
    if (v == 1)
    {
        ds >> sym.stringrep;
    
        //get the ID number for this symbol
        sym.id = Symbol::getNewID(sym.stringrep);
    }
    else
        throw version_error(v, "1", r_symbol, CODELOC);
    
    return ds;
}

/** Null constructor */
Symbol::Symbol() : ExBase(), id(0), stringrep(QString::null)
{}

/** Construct a new symbol, with string representation 'rep' */
Symbol::Symbol(const QString &rep) 
       : ExBase(), id( Symbol::getNewID(rep) ), stringrep(rep)
{}

/** Copy constructor */
Symbol::Symbol(const Symbol &other)
       : ExBase(), id(other.id), stringrep(other.stringrep)
{}

/** Destructor */
Symbol::~Symbol()
{}

/** Assignment operator */
Symbol& Symbol::operator=(const Symbol &other)
{
    id = other.id;
    stringrep = other.stringrep;
    return *this;
}

/** Comparison operator */
bool Symbol::operator==(const ExBase &other) const
{
    const Symbol *sym = dynamic_cast<const Symbol*>(&other);
    
    return sym != 0 and typeid(other).name() == typeid(*this).name()
            and sym->ID() == this->ID();
}

/** Return a string representation of this symbol */
QString Symbol::toString() const
{
    return stringrep;
}

/** There are no child expressions in a symbol */
Expressions Symbol::children() const
{
    return Expressions();
}

/** Evaluate this symbol - returns the value of the symbol in 'values' if
    it is present, else it returns 0.0 */
double Symbol::evaluate(const Values &values) const
{
    return values.value(*this);
}

/** Evaluate this symbol - returns the value of the symbol in 'values' if
    it is present, else it returns 0 */
Complex Symbol::evaluate(const ComplexValues &values) const
{
    return values.value(*this);
}

/** Differentiate this symbol with respect to 'sym'. This returns 1.0 if this
    is 'sym', else it returns 0.0 */
Expression Symbol::differentiate(const Symbol &sym) const
{
    //assume int(true) == 1 and int(false) == 0
    return Expression( int( sym.ID() == ID() ) );
}

/** Integrate this symbol with respect to 'sym'. If 'sym' == this, then 
    return 0.5 sym^2, else return *this * sym */
Expression Symbol::integrate(const Symbol &sym) const
{
    if (sym.ID() == ID())
        return 0.5 * pow(sym,2);
    else
        return *this * sym;
}

/** Return the expression that matches this symbol in 'identities' - or return  
    an expression holding only this symbol if it does no exist in 'identities' */
Expression Symbol::substitute(const Identities &identities) const
{
    return identities.expression(*this);
}

/** Return this symbol */
Symbols Symbol::symbols() const
{
    return Symbols(*this);
}

/** This is not a function */
Functions Symbol::functions() const
{
    return Functions();
}

/** Is this a function of 'symbol' */
bool Symbol::isFunction(const Symbol &sym) const
{
    return (id != 0) and (sym.ID() == ID());
}

/** A symbol is by definition not constant */
bool Symbol::isConstant() const
{
    return false;
}
