#include <cmath>

#include <boost/variant.hpp>

#include "src/storage/expressions/UnaryNumericalFunctionExpression.h"
#include "src/storage/expressions/IntegerLiteralExpression.h"
#include "src/storage/expressions/DoubleLiteralExpression.h"
#include "src/utility/macros.h"
#include "src/exceptions/InvalidTypeException.h"

namespace storm {
    namespace expressions {
        UnaryNumericalFunctionExpression::UnaryNumericalFunctionExpression(ExpressionManager const& manager, Type const& type, std::shared_ptr<BaseExpression const> const& operand, OperatorType operatorType) : UnaryExpression(manager, type, operand), operatorType(operatorType) {
            // Intentionally left empty.
        }
        
        UnaryNumericalFunctionExpression::OperatorType UnaryNumericalFunctionExpression::getOperatorType() const {
            return this->operatorType;
        }
        
        storm::expressions::OperatorType UnaryNumericalFunctionExpression::getOperator() const {
            storm::expressions::OperatorType result = storm::expressions::OperatorType::Minus;
            switch (this->getOperatorType()) {
                case OperatorType::Minus: result = storm::expressions::OperatorType::Minus; break;
                case OperatorType::Floor: result = storm::expressions::OperatorType::Floor; break;
                case OperatorType::Ceil: result = storm::expressions::OperatorType::Ceil; break;
            }
            return result;
        }
        
        int_fast64_t UnaryNumericalFunctionExpression::evaluateAsInt(Valuation const* valuation) const {
            STORM_LOG_THROW(this->hasIntegerType(), storm::exceptions::InvalidTypeException, "Unable to evaluate expression as integer.");

            int_fast64_t result = this->getOperand()->evaluateAsInt(valuation);
            switch (this->getOperatorType()) {
                case OperatorType::Minus: result = -result; break;
                case OperatorType::Floor: result = std::floor(result); break;
                case OperatorType::Ceil: result = std::ceil(result); break;
            }
            return result;
        }
        
        double UnaryNumericalFunctionExpression::evaluateAsDouble(Valuation const* valuation) const {
            STORM_LOG_THROW(this->hasNumericalType(), storm::exceptions::InvalidTypeException, "Unable to evaluate expression as double.");

            double result = this->getOperand()->evaluateAsDouble(valuation);
            switch (this->getOperatorType()) {
                case OperatorType::Minus: result = -result; break;
                case OperatorType::Floor: result = std::floor(result); break;
                case OperatorType::Ceil: result = std::ceil(result); break;
            }
            return result;
        }
        
        std::shared_ptr<BaseExpression const> UnaryNumericalFunctionExpression::simplify() const {
            std::shared_ptr<BaseExpression const> operandSimplified = this->getOperand()->simplify();
            
            if (operandSimplified->isLiteral()) {
                boost::variant<int_fast64_t, double> operandEvaluation;

                if (operandSimplified->hasIntegerType()) {
                    operandEvaluation = operandSimplified->evaluateAsInt();
                } else {
                    operandEvaluation = operandSimplified->evaluateAsDouble();
                }
                
                boost::variant<int_fast64_t, double> result;
                if (operandSimplified->hasIntegerType()) {
                    int_fast64_t value = 0;
                    switch (this->getOperatorType()) {
                        case OperatorType::Minus: value = -boost::get<int_fast64_t>(operandEvaluation); break;
                        case OperatorType::Floor: value = std::floor(boost::get<int_fast64_t>(operandEvaluation)); break;
                        case OperatorType::Ceil: value = std::ceil(boost::get<int_fast64_t>(operandEvaluation)); break;
                    }
                    return std::shared_ptr<BaseExpression>(new IntegerLiteralExpression(this->getManager(), value));
                } else {
                    double value = 0;
                    switch (this->getOperatorType()) {
                        case OperatorType::Minus: value = -boost::get<double>(operandEvaluation); break;
                        case OperatorType::Floor: value = std::floor(boost::get<double>(operandEvaluation)); break;
                        case OperatorType::Ceil: value = std::ceil(boost::get<double>(operandEvaluation)); break;
                    }
                    return std::shared_ptr<BaseExpression>(new DoubleLiteralExpression(this->getManager(), value));
                }
            }
            
            if (operandSimplified.get() == this->getOperand().get()) {
                return this->shared_from_this();
            } else {
                return std::shared_ptr<BaseExpression>(new UnaryNumericalFunctionExpression(this->getManager(), this->getType(), operandSimplified, this->getOperatorType()));
            }
        }
        
        boost::any UnaryNumericalFunctionExpression::accept(ExpressionVisitor& visitor) const {
            return visitor.visit(*this);
        }
        
        void UnaryNumericalFunctionExpression::printToStream(std::ostream& stream) const {
            switch (this->getOperatorType()) {
                case OperatorType::Minus: stream << "-("; break;
                case OperatorType::Floor: stream << "floor("; break;
                case OperatorType::Ceil: stream << "ceil("; break;
            }
            stream << *this->getOperand() << ")";
        }
    }
}