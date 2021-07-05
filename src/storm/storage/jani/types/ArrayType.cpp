#include "storm/storage/jani/types/ArrayType.h"

namespace storm {
    namespace jani {
        ArrayType::ArrayType(JaniType* childType) : JaniType(), childType(childType){
            // Intentionally left empty
        }

        bool ArrayType::isArrayType() const {
            return true;
        }

        bool ArrayType::isBoundedType() const {
            return childType->isBoundedType();
        }

        JaniType* ArrayType::getChildType() const {
            return childType;
        }

        std::string ArrayType::getStringRepresentation() const {
            return "array[" + getChildType()->getStringRepresentation() + "]";
        }

        void ArrayType::setLowerBound(storm::expressions::Expression const& expression) {
            childType->setLowerBound(expression);
        }

        void ArrayType::setUpperBound(storm::expressions::Expression const& expression) {
            childType->setUpperBound(expression);
        }

        storm::expressions::Expression const& ArrayType::getLowerBound() const {
            return childType->getLowerBound();
        }

        storm::expressions::Expression const& ArrayType::getUpperBound() const {
            return childType->getUpperBound();
        }
    }
}