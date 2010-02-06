// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Node.pypp.hpp"

namespace bp = boost::python;

#include "Siren/errors.h"

#include "node.h"

#include "nodes.h"

#include "promise.h"

#include "promises.h"

#include "workpacket.h"

#include "workqueue.h"

#include "node.h"

#include "Siren/logger.h"

SireCluster::Node __copy__(const SireCluster::Node &other){ return SireCluster::Node(other); }

#include "Siren/str.hpp"

void register_Node_class(){

    { //::SireCluster::Node
        typedef bp::class_< SireCluster::Node, bp::bases< Siren::Handle > > Node_exposer_t;
        Node_exposer_t Node_exposer = Node_exposer_t( "Node" );
        bp::scope Node_scope( Node_exposer );
        Node_exposer.def( bp::init< >() );
        Node_exposer.def( bp::init< SireCluster::Node const & >(( bp::arg("other") )) );
        { //::SireCluster::Node::hashCode
        
            typedef ::uint ( ::SireCluster::Node::*hashCode_function_type )(  ) const;
            hashCode_function_type hashCode_function_value( &::SireCluster::Node::hashCode );
            
            Node_exposer.def( 
                "hashCode"
                , hashCode_function_value );
        
        }
        { //::SireCluster::Node::isBusy
        
            typedef bool ( ::SireCluster::Node::*isBusy_function_type )(  ) const;
            isBusy_function_type isBusy_function_value( &::SireCluster::Node::isBusy );
            
            Node_exposer.def( 
                "isBusy"
                , isBusy_function_value );
        
        }
        { //::SireCluster::Node::isFree
        
            typedef bool ( ::SireCluster::Node::*isFree_function_type )(  ) const;
            isFree_function_type isFree_function_value( &::SireCluster::Node::isFree );
            
            Node_exposer.def( 
                "isFree"
                , isFree_function_value );
        
        }
        { //::SireCluster::Node::merge
        
            typedef ::SireCluster::Nodes ( *merge_function_type )( ::SireCluster::Nodes,::SireCluster::Nodes );
            merge_function_type merge_function_value( &::SireCluster::Node::merge );
            
            Node_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("nodes0"), bp::arg("nodes1") ) );
        
        }
        { //::SireCluster::Node::merge
        
            typedef ::SireCluster::Nodes ( *merge_function_type )( ::SireCluster::Node,::SireCluster::Node );
            merge_function_type merge_function_value( &::SireCluster::Node::merge );
            
            Node_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("node0"), bp::arg("node1") ) );
        
        }
        { //::SireCluster::Node::merge
        
            typedef ::SireCluster::Nodes ( *merge_function_type )( ::SireCluster::Nodes,::SireCluster::Node );
            merge_function_type merge_function_value( &::SireCluster::Node::merge );
            
            Node_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("nodes0"), bp::arg("node1") ) );
        
        }
        { //::SireCluster::Node::merge
        
            typedef ::SireCluster::Nodes ( *merge_function_type )( ::SireCluster::Node,::SireCluster::Nodes );
            merge_function_type merge_function_value( &::SireCluster::Node::merge );
            
            Node_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("node0"), bp::arg("nodes1") ) );
        
        }
        Node_exposer.def( bp::self != bp::self );
        { //::SireCluster::Node::operator=
        
            typedef ::SireCluster::Node & ( ::SireCluster::Node::*assign_function_type )( ::SireCluster::Node const & ) ;
            assign_function_type assign_function_value( &::SireCluster::Node::operator= );
            
            Node_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        Node_exposer.def( bp::self == bp::self );
        { //::SireCluster::Node::submit
        
            typedef ::SireCluster::Promise ( ::SireCluster::Node::*submit_function_type )( ::SireCluster::WorkPacket const & ) ;
            submit_function_type submit_function_value( &::SireCluster::Node::submit );
            
            Node_exposer.def( 
                "submit"
                , submit_function_value
                , ( bp::arg("workpacket") ) );
        
        }
        { //::SireCluster::Node::submit
        
            typedef ::SireCluster::Promises ( ::SireCluster::Node::*submit_function_type )( ::QList< Siren::ObjPtr< SireCluster::WorkPacket > > const & ) ;
            submit_function_type submit_function_value( &::SireCluster::Node::submit );
            
            Node_exposer.def( 
                "submit"
                , submit_function_value
                , ( bp::arg("workpackets") ) );
        
        }
        { //::SireCluster::Node::submit
        
            typedef ::SireCluster::Promises ( ::SireCluster::Node::*submit_function_type )( ::QVector< Siren::ObjPtr< SireCluster::WorkPacket > > const & ) ;
            submit_function_type submit_function_value( &::SireCluster::Node::submit );
            
            Node_exposer.def( 
                "submit"
                , submit_function_value
                , ( bp::arg("workpackets") ) );
        
        }
        { //::SireCluster::Node::toString
        
            typedef ::QString ( ::SireCluster::Node::*toString_function_type )(  ) const;
            toString_function_type toString_function_value( &::SireCluster::Node::toString );
            
            Node_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        Node_exposer.staticmethod( "merge" );
        Node_exposer.def( "__copy__", &__copy__);
        Node_exposer.def( "__deepcopy__", &__copy__);
        Node_exposer.def( "clone", &__copy__);
        Node_exposer.def( "__str__", &__str__< ::SireCluster::Node > );
        Node_exposer.def( "__repr__", &__str__< ::SireCluster::Node > );
    }

}