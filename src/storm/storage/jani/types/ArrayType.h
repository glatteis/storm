#pragma once

#include "JaniType.h"

namespace storm {
    namespace jani {
        class ArrayType : public JaniType {
        public:
            ArrayType(JaniType* childType);
            bool isArrayType() const override;
            bool isBoundedType() const override;
            JaniType* getChildType() const override;
            std::string getStringRepresentation() const override;
            void setLowerBound(storm::expressions::Expression const& expression) override;
            void setUpperBound(storm::expressions::Expression const& expression) override;
            storm::expressions::Expression const& getLowerBound() const override;
            storm::expressions::Expression const& getUpperBound() const override;

        private:
            JaniType* childType;

        };
    }
}

