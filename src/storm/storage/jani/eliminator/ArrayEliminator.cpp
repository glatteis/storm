#include "ArrayEliminator.h"

#include <unordered_map>

#include "storm/storage/expressions/ExpressionVisitor.h"
#include "storm/storage/jani/visitor/JaniExpressionVisitor.h"
#include "storm/storage/jani/Variable.h"
#include "storm/storage/jani/Model.h"
#include "storm/storage/jani/Property.h"
#include "storm/storage/jani/traverser/JaniTraverser.h"
#include "storm/storage/jani/traverser/ArrayExpressionFinder.h"

#include "storm/storage/expressions/Expressions.h"
#include "storm/storage/jani/expressions/JaniExpressions.h"
#include "storm/storage/expressions/ExpressionManager.h"

#include "storm/exceptions/NotSupportedException.h"
#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/OutOfRangeException.h"

namespace storm {
    
    
    
    namespace jani {
        namespace detail {
            
            class MaxArraySizeExpressionVisitor : public storm::expressions::ExpressionVisitor, public storm::expressions::JaniExpressionVisitor {
            public:
                MaxArraySizeExpressionVisitor() = default;
                virtual ~MaxArraySizeExpressionVisitor() = default;

                std::size_t getMaxSize(storm::expressions::Expression const& expression, std::unordered_map<storm::expressions::Variable, std::size_t> const& arrayVariableSizeMap) {
                    return boost::any_cast<std::size_t>(expression.accept(*this, &arrayVariableSizeMap));
                }

                std::size_t getMaxSizeAt(storm::expressions::Expression const& expression, std::unordered_map<storm::expressions::Variable, std::size_t> const& arrayVariableSizeMap, int const index) {
                    std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int> dataPair = {arrayVariableSizeMap, index};
                    return boost::any_cast<std::size_t>(expression.accept(*this, dataPair));
                }
     
                virtual boost::any visit(storm::expressions::IfThenElseExpression const& expression, boost::any const& data) override {
                    if (expression.getCondition()->containsVariables()) {
                        return std::max<std::size_t>(boost::any_cast<std::size_t>(expression.getThenExpression()->accept(*this, data)), boost::any_cast<std::size_t>(expression.getElseExpression()->accept(*this, data)));
                    } else {
                        if (expression.getCondition()->evaluateAsBool()) {
                            return boost::any_cast<std::size_t>(expression.getThenExpression()->accept(*this, data));
                        }
                        return boost::any_cast<std::size_t>(expression.getElseExpression()->accept(*this, data));
                    }
                }
                
                virtual boost::any visit(storm::expressions::BinaryBooleanFunctionExpression const& expression, boost::any const& data) override {
                    return std::max<std::size_t>(boost::any_cast<std::size_t>(expression.getFirstOperand()->accept(*this, data)), boost::any_cast<std::size_t>(expression.getSecondOperand()->accept(*this, data)));
                }
                
                virtual boost::any visit(storm::expressions::BinaryNumericalFunctionExpression const& expression, boost::any const& data) override {
                    return std::max<std::size_t>(boost::any_cast<std::size_t>(expression.getFirstOperand()->accept(*this, data)), boost::any_cast<std::size_t>(expression.getSecondOperand()->accept(*this, data)));
                }
                
                virtual boost::any visit(storm::expressions::BinaryRelationExpression const& expression, boost::any const& data) override {
                    return std::max<std::size_t>(boost::any_cast<std::size_t>(expression.getFirstOperand()->accept(*this, data)), boost::any_cast<std::size_t>(expression.getSecondOperand()->accept(*this, data)));
                }
                
                virtual boost::any visit(storm::expressions::VariableExpression const& expression, boost::any const& data) override {
                    if (data.type() == typeid(std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>)) {
                        auto arrayVariableSizeMap = boost::any_cast<std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>>(data);
                        if (expression.getType().isArrayType()) {
                            // varIt is entry from map std::unordered_map<storm::expressions::Variable
                            auto varIt = arrayVariableSizeMap.first.find(expression.getVariable());
                            if (varIt != arrayVariableSizeMap.first.end()) {
                                // arrayVariableSize.Map.second is the index
                                return varIt->first.getArraySize(arrayVariableSizeMap.second);
                            }
                        }
                    } else if (data.type() == typeid(std::unordered_map<storm::expressions::Variable, std::size_t>)) {
                        auto arrayVariableSizeMap = boost::any_cast<std::unordered_map<storm::expressions::Variable, std::size_t>>(data);
                        if (expression.getType().isArrayType()) {
                            auto varIt = arrayVariableSizeMap.find(expression.getVariable());
                            if (varIt != arrayVariableSizeMap.end()) {
                                return varIt->second;
                            }
                        }
                    }
                    return static_cast<std::size_t>(0);
                }
                
                virtual boost::any visit(storm::expressions::UnaryBooleanFunctionExpression const& expression, boost::any const& data) override {
                    return boost::any_cast<std::size_t>(expression.getOperand()->accept(*this, data));
                }
                
                virtual boost::any visit(storm::expressions::UnaryNumericalFunctionExpression const& expression, boost::any const& data) override {
                    return boost::any_cast<std::size_t>(expression.getOperand()->accept(*this, data));
                }
                
                virtual boost::any visit(storm::expressions::BooleanLiteralExpression const&, boost::any const&) override {
                    return 0;
                }
                
                virtual boost::any visit(storm::expressions::IntegerLiteralExpression const&, boost::any const&) override {
                    return 0;
                }
                
                virtual boost::any visit(storm::expressions::RationalLiteralExpression const&, boost::any const&) override {
                    return 0;
                }
                
                virtual boost::any visit(storm::expressions::ValueArrayExpression const& expression, boost::any const& data) override {
                    if (data.type() == typeid(std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>)) {
                        // This is to get the size of the array at nested array number i
                        auto dataPair = boost::any_cast<std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>>(data);
                        return visit(expression.getElements(), dataPair.second);
                    } else {
                        STORM_LOG_ASSERT(expression.size()->isIntegerLiteralExpression(), "unexpected kind of size expression of ValueArrayExpression (" << expression.size()->toExpression() << ").");
                    }
                    return static_cast<std::size_t>(expression.size()->evaluateAsInt());
                }

                virtual boost::any visit(storm::expressions::ValueArrayExpression::ValueArrayElements const& expression, boost::any const& data) override {
                    if (data.type() == typeid(int)) {
                        // This is to get the size of the array at nested array number i
                        auto i = boost::any_cast<int>(data);
                        return expression.getSizes().at(i);
                    } else {
                        assert (false);
                        // TODO: implement
                    }
                    return static_cast<std::size_t>(0);
                }
                
                virtual boost::any visit(storm::expressions::ConstructorArrayExpression const& expression, boost::any const& data) override {
                    if (data.type() == typeid(std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>)) {
                        auto dataPair = boost::any_cast<std::pair<std::unordered_map<storm::expressions::Variable, std::size_t>, int>>(data);
                        return static_cast<std::size_t>(expression.size(dataPair.second)->evaluateAsInt());
                    } else {
                        if (!expression.size()->containsVariables()) {
                            return static_cast<std::size_t>(expression.size()->evaluateAsInt());
                        } else {
                            auto vars = expression.size()->toExpression().getVariables();
                            std::string variables = "";
                            for (auto const &v : vars) {
                                if (variables != "") {
                                    variables += ", ";
                                }
                                variables += v.getName();
                            }
                            if (vars.size() == 1) {
                                variables = "variable " + variables;
                            } else {
                                variables = "variables " + variables;
                            }
                            STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unable to determine array size: Size of ConstructorArrayExpression '" << expression << "' still contains the " << variables << ".");
                            return static_cast<std::size_t>(0);
                        }
                    }
                }
                
                virtual boost::any visit(storm::expressions::ArrayAccessExpression const& expression, boost::any const& data) override {
                    return expression.getSecondOperand()->accept(*this, data);
                }

                virtual boost::any visit(storm::expressions::ArrayAccessIndexExpression const& expression, boost::any const& data) override {
                    if (expression.getFirstOperand() == expression.getSecondOperand()) {
                        return expression.getFirstOperand()->accept(*this, data);
                    } else {
                        auto firstPart = expression.getFirstOperand()->accept(*this, data);
                        assert (false);
                    }
                }

                virtual boost::any visit(storm::expressions::FunctionCallExpression const&, boost::any const&) override {
                    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Found Function call expression within an array expression. This is not expected since functions are expected to be eliminated at this point.");
                    return 0;
                }
                
            };
            
            class ArrayExpressionEliminationVisitor : public storm::expressions::ExpressionVisitor, public storm::expressions::JaniExpressionVisitor {
            public:
                
                typedef std::shared_ptr<storm::expressions::BaseExpression const> BaseExprPtr;
                class ResultType {
                public:
                    ResultType(ResultType const& other) = default;
                    ResultType(BaseExprPtr expression) : expression(expression), arrayOutOfBoundsMessage("") {}
                    ResultType(std::string arrayOutOfBoundsMessage) : expression(nullptr), arrayOutOfBoundsMessage(arrayOutOfBoundsMessage) {}
                    BaseExprPtr& expr() {
                        STORM_LOG_ASSERT(!isArrayOutOfBounds(), "Tried to get the result expression, but " << arrayOutOfBoundsMessage);
                        return expression;
                    };
                    bool isArrayOutOfBounds() { return arrayOutOfBoundsMessage != ""; };
                    std::string const& outOfBoundsMessage() const { return arrayOutOfBoundsMessage; }
                private:
                    BaseExprPtr expression;
                    std::string arrayOutOfBoundsMessage;
                };
                
                ArrayExpressionEliminationVisitor(std::unordered_map<storm::expressions::Variable, std::vector<storm::jani::Variable const*>> const& replacements, std::unordered_map<storm::expressions::Variable, std::size_t> const& sizes) : replacements(replacements), arraySizes(sizes) {}
                virtual ~ArrayExpressionEliminationVisitor() = default;
    
                std::vector<storm::expressions::Expression> eliminate(std::vector<storm::expressions::Expression> const& expressions) {
                    // here, data is the accessed index of the most recent array access expression. Initially, there is none.
                    // TODO: implement this for the expressions, these are index expressions, need to be mapped to 1 index
                    std::vector<storm::expressions::Expression> result;
                    for (auto& expression : expressions) {
                        result.push_back(eliminate(expression));
                        assert(!containsArrayExpression(result.back()));
                    }
                    return result;
                }

                storm::expressions::Expression eliminate(storm::expressions::Expression const& expression) {
                    // here, data is the accessed index of the most recent array access expression. Initially, there is none.
                    auto res = boost::any_cast<ResultType>(expression.accept(*this, boost::any()));
                    STORM_LOG_THROW(!res.isArrayOutOfBounds(), storm::exceptions::OutOfRangeException, res.outOfBoundsMessage());
                    STORM_LOG_ASSERT(!containsArrayExpression(res.expr()->toExpression()), "Expression still contains array expressions. Before: " << std::endl << expression << std::endl << "After:" << std::endl << res.expr()->toExpression());
                    return res.expr()->simplify();
                }
                
                virtual boost::any visit(storm::expressions::IfThenElseExpression const& expression, boost::any const& data) override {
                    // for the condition expression, outer array accesses should not matter.
                    ResultType conditionResult = boost::any_cast<ResultType>(expression.getCondition()->accept(*this, boost::any()));
                    if (conditionResult.isArrayOutOfBounds()) {
                        return conditionResult;
                    }
                    
                    // We need to handle expressions of the kind '42<size : A[42] : 0', where size is a variable and A[42] might be out of bounds.
                    ResultType thenResult = boost::any_cast<ResultType>(expression.getThenExpression()->accept(*this, data));
                    ResultType elseResult = boost::any_cast<ResultType>(expression.getElseExpression()->accept(*this, data));
                    if (thenResult.isArrayOutOfBounds()) {
                        if (elseResult.isArrayOutOfBounds()) {
                            return ResultType(thenResult.outOfBoundsMessage() + " and " + elseResult.outOfBoundsMessage());
                        } else {
                            // Assume the else expression
                            return elseResult;
                        }
                    } else if (elseResult.isArrayOutOfBounds()) {
                        // Assume the then expression
                        return thenResult;
                    } else {
                        // If the arguments did not change, we simply push the expression itself.
                        if (conditionResult.expr().get() == expression.getCondition().get() && thenResult.expr().get() == expression.getThenExpression().get() && elseResult.expr().get() == expression.getElseExpression().get()) {
                            return ResultType(expression.getSharedPointer());
                        } else {
                            return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::IfThenElseExpression(expression.getManager(), thenResult.expr()->getType(), conditionResult.expr(), thenResult.expr(), elseResult.expr()))));
                        }
                    }
                }
        
                virtual boost::any visit(storm::expressions::BinaryBooleanFunctionExpression const& expression, boost::any const& data) override {
                    STORM_LOG_ASSERT(data.empty(), "BinaryBooleanFunctionExpressions should not be direct subexpressions of array access expressions. However, the expression " << expression << " is.");
                    ResultType firstResult = boost::any_cast<ResultType>(expression.getFirstOperand()->accept(*this, data));
                    ResultType secondResult = boost::any_cast<ResultType>(expression.getSecondOperand()->accept(*this, data));
                    
                    if (firstResult.isArrayOutOfBounds()) {
                        return firstResult;
                    } else if (secondResult.isArrayOutOfBounds()) {
                        return secondResult;
                    }
        
                    // If the arguments did not change, we simply push the expression itself.
                    if (firstResult.expr().get() == expression.getFirstOperand().get() && secondResult.expr().get() == expression.getSecondOperand().get()) {
                        return ResultType(expression.getSharedPointer());
                    } else {
                        return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::BinaryBooleanFunctionExpression(expression.getManager(), expression.getType(), firstResult.expr(), secondResult.expr(), expression.getOperatorType()))));
                    }
                }
        
                virtual boost::any visit(storm::expressions::BinaryNumericalFunctionExpression const& expression, boost::any const& data) override {
                    STORM_LOG_ASSERT(data.empty(), "BinaryNumericalFunctionExpression should not be direct subexpressions of array access expressions. However, the expression " << expression << " is.");
                    ResultType firstResult = boost::any_cast<ResultType>(expression.getFirstOperand()->accept(*this, data));
                    ResultType secondResult = boost::any_cast<ResultType>(expression.getSecondOperand()->accept(*this, data));
                    
                    if (firstResult.isArrayOutOfBounds()) {
                        return firstResult;
                    } else if (secondResult.isArrayOutOfBounds()) {
                        return secondResult;
                    }
                    
                    // If the arguments did not change, we simply push the expression itself.
                    if (firstResult.expr().get() == expression.getFirstOperand().get() && secondResult.expr().get() == expression.getSecondOperand().get()) {
                        return ResultType(expression.getSharedPointer());
                    } else {
                        return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::BinaryNumericalFunctionExpression(expression.getManager(), expression.getType(), firstResult.expr(), secondResult.expr(), expression.getOperatorType()))));
                    }
                }
        
                virtual boost::any visit(storm::expressions::BinaryRelationExpression const& expression, boost::any const& data) override {
                    STORM_LOG_ASSERT(data.empty(), "BinaryRelationExpression should not be direct subexpressions of array access expressions. However, the expression " << expression << " is.");
                    ResultType firstResult = boost::any_cast<ResultType>(expression.getFirstOperand()->accept(*this, data));
                    ResultType secondResult = boost::any_cast<ResultType>(expression.getSecondOperand()->accept(*this, data));
                    
                    if (firstResult.isArrayOutOfBounds()) {
                        return firstResult;
                    } else if (secondResult.isArrayOutOfBounds()) {
                        return secondResult;
                    }
                    
                    // If the arguments did not change, we simply push the expression itself.
                    if (firstResult.expr().get() == expression.getFirstOperand().get() && secondResult.expr().get() == expression.getSecondOperand().get()) {
                        return ResultType(expression.getSharedPointer());
                    } else {
                        return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::BinaryRelationExpression(expression.getManager(), expression.getType(), firstResult.expr(), secondResult.expr(), expression.getRelationType()))));
                    }
                }
                
                virtual boost::any visit(storm::expressions::VariableExpression const& expression, boost::any const& data) override {
                    if (expression.getType().isArrayType()) {
                        STORM_LOG_THROW(!data.empty(), storm::exceptions::NotSupportedException, "Unable to translate array variable to basic variable, since it does not seem to be within an array access expression.");
                        uint64_t index;
                        if (data.type() == typeid(uint64_t)) {
                            index = boost::any_cast<uint64_t>(data);
                        } else {
                            assert (data.type() == typeid(ResultType));
                            auto resultType = boost::any_cast<ResultType>(data);
                            if (resultType.expr()->containsVariables()) {
                                assert (false);
                            } else {
                                index = resultType.expr()->evaluateAsInt();
                            }
                        }
                        STORM_LOG_ASSERT(replacements.find(expression.getVariable()) != replacements.end(), "Unable to find array variable " << expression << " in array replacements.");
                        auto const& arrayVarReplacements = replacements.at(expression.getVariable());
                        if (index >= arrayVarReplacements.size()) {
                            return ResultType("Array index " + std::to_string(index) + " for variable " + expression.getVariableName() + " is out of bounds.");
                        }
                        return ResultType(arrayVarReplacements[index]->getExpressionVariable().getExpression().getBaseExpressionPointer());
                    } else {
                        STORM_LOG_ASSERT(data.empty(), "VariableExpression of non-array variable should not be a subexpressions of array access expressions. However, the expression " << expression << " is.");
                        return ResultType(expression.getSharedPointer());
                    }
                }
        
                virtual boost::any visit(storm::expressions::UnaryBooleanFunctionExpression const& expression, boost::any const& data) override {
                    STORM_LOG_ASSERT(data.empty(), "UnaryBooleanFunctionExpression should not be direct subexpressions of array access expressions. However, the expression " << expression << " is.");
                    ResultType operandResult = boost::any_cast<ResultType>(expression.getOperand()->accept(*this, data));
                    
                    if (operandResult.isArrayOutOfBounds()) {
                        return operandResult;
                    }
                    
                    // If the argument did not change, we simply push the expression itself.
                    if (operandResult.expr().get() == expression.getOperand().get()) {
                        return ResultType(expression.getSharedPointer());
                    } else {
                        return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::UnaryBooleanFunctionExpression(expression.getManager(), expression.getType(), operandResult.expr(), expression.getOperatorType()))));
                    }
                }
        
                virtual boost::any visit(storm::expressions::UnaryNumericalFunctionExpression const& expression, boost::any const& data) override {
                    STORM_LOG_ASSERT(data.empty(), "UnaryBooleanFunctionExpression should not be direct subexpressions of array access expressions. However, the expression " << expression << " is.");
        
                    ResultType operandResult = boost::any_cast<ResultType>(expression.getOperand()->accept(*this, data));
                    
                    if (operandResult.isArrayOutOfBounds()) {
                        return operandResult;
                    }
                   
                    // If the argument did not change, we simply push the expression itself.
                    if (operandResult.expr().get() == expression.getOperand().get()) {
                        return ResultType(expression.getSharedPointer());
                    } else {
                        return ResultType(std::const_pointer_cast<storm::expressions::BaseExpression const>(std::shared_ptr<storm::expressions::BaseExpression>(new storm::expressions::UnaryNumericalFunctionExpression(expression.getManager(), expression.getType(), operandResult.expr(), expression.getOperatorType()))));
                    }
                }
                
                virtual boost::any visit(storm::expressions::BooleanLiteralExpression const& expression, boost::any const&) override {
                    return ResultType(expression.getSharedPointer());
                }
                
                virtual boost::any visit(storm::expressions::IntegerLiteralExpression const& expression, boost::any const&) override {
                    return ResultType(expression.getSharedPointer());
                }
                
                virtual boost::any visit(storm::expressions::RationalLiteralExpression const& expression, boost::any const&) override {
                    return ResultType(expression.getSharedPointer());
                }
                
                virtual boost::any visit(storm::expressions::ValueArrayExpression const& expression, boost::any const& data) override {
                    STORM_LOG_THROW(!data.empty(), storm::exceptions::NotSupportedException, "Unable to translate ValueArrayExpression to element expression since it does not seem to be within an array access expression.");
                    uint64_t index = boost::any_cast<uint64_t>(data);
                    STORM_LOG_ASSERT(expression.size()->isIntegerLiteralExpression(), "unexpected kind of size expression of ValueArrayExpression (" << expression.size()->toExpression() << ").");
                    if (index >= static_cast<uint64_t>(expression.size()->evaluateAsInt())) {
                        return ResultType("Array index " + std::to_string(index) + " for ValueArrayExpression " + expression.toExpression().toString() + " is out of bounds.");
                    }
                    return ResultType(boost::any_cast<ResultType>(expression.at(index)->accept(*this, boost::any())));
                }

                virtual boost::any visit(storm::expressions::ValueArrayExpression::ValueArrayElements const& expression, boost::any const& data) override {
                    assert (false); // Left empty as we should not visit this one, we handle this with expression.at in the visit for ValueArrayExpression
                }

                virtual boost::any visit(storm::expressions::ConstructorArrayExpression const& expression, boost::any const& data) override {
                    STORM_LOG_THROW(!data.empty(), storm::exceptions::NotSupportedException, "Unable to translate ValueArrayExpression to element expression since it does not seem to be within an array access expression.");
                    uint64_t index;
                    if (data.type() == typeid(uint64_t)) {
                        index = boost::any_cast<uint64_t>(data);
                    } else {
                        auto resultType = boost::any_cast<ResultType>(data);
                        if (resultType.expr()->containsVariables()) {
                            assert (false);
                        } else {
                            index = resultType.expr()->evaluateAsInt();
                        }
                    }
                    if (expression.size()->containsVariables()) {
                        STORM_LOG_WARN("Ignoring length of constructorArrayExpression " << expression << " as it still contains variables.");
                    } else {
                        if (index >= static_cast<uint64_t>(expression.size()->evaluateAsInt())) {
                            return ResultType("Array index " + std::to_string(index) + " for ConstructorArrayExpression " + expression.toExpression().toString() + " is out of bounds.");
                        }
                    }
                    return ResultType(boost::any_cast<ResultType>(expression.at(index)->accept(*this, boost::any())));
                }
                
                virtual boost::any visit(storm::expressions::ArrayAccessExpression const& expression, boost::any const&) override {
                    // We let the data consist of the pointer referring to the arrayAccessExpression variable
                    // and an int telling us the arrayNumber of the current arrayAccessIndex
                    // e.g. for a[1][4] the data will be a pointer to a and the arrayNumber will be 0, telling us we are considering the first part of the array
                    std::pair<std::shared_ptr<storm::expressions::BaseExpression const>, uint_fast64_t> data = {expression.getFirstOperand(), 0};
                    return expression.getSecondOperand()->accept(*this, data);
                }

                virtual boost::any visit(storm::expressions::ArrayAccessIndexExpression const& expression, boost::any const& data) override {
                    // e.g. for a[1][4] the data is a pointer to a
                    // if the arrayNumber 0 we are considering 1 in a[1]
                    // if the arrayNumber is 1 we are considering the 4 in a[4]
                    assert (data.type() == typeid(std::pair<std::shared_ptr<storm::expressions::BaseExpression const>, uint_fast64_t>));

                    auto castedData = boost::any_cast<std::pair<std::shared_ptr<storm::expressions::BaseExpression const>, uint_fast64_t>>(data);
                    auto arrayExpression = castedData.first;
                    auto arrayNumber = castedData.second;
                    // Go over other part of array, the first part of resultSecond will be a map containing of index and isCurrentIndex expressions, the second part is the size of the array so far
                    // e.g. array(4, array(5, array(3))) will give us for:
                    // - the array(3): 0
                    // - array(5, array(3)): 3
                    // - array(4, array(5, array(3))): 15
                    std::pair<std::shared_ptr<storm::expressions::BaseExpression const>, uint_fast64_t> newData = {arrayExpression, arrayNumber + 1};
                    boost::optional<std::pair<std::map<uint_fast64_t, storm::expressions::Expression>, uint_fast64_t>> resultSecond;
                    auto sizeSoFar = 1;
                    if (expression.getFirstOperand() != expression.getSecondOperand()) {
                        resultSecond = boost::any_cast<std::pair<std::map<uint_fast64_t, storm::expressions::Expression>, uint_fast64_t> >(expression.getSecondOperand()->accept(*this, newData));
                        sizeSoFar = resultSecond->second;
                    }

                    uint64_t size = MaxArraySizeExpressionVisitor().getMaxSizeAt(arrayExpression->toExpression(), arraySizes, arrayNumber);
                    STORM_LOG_THROW(size > 0, storm::exceptions::NotSupportedException, "Unable to get size of array expression for array access " << expression << ".");

                    std::map<uint_fast64_t, storm::expressions::Expression> expressions;

                    if (expression.getFirstOperand() == expression.getSecondOperand()) {
                        // Last array, so resultSecond will not be initialized
                        if (expression.getFirstOperand()->containsVariables()) {
                            for (size_t index = 0; index < size; ++index) {
                                storm::expressions::Expression isCurrentIndex = boost::any_cast<ResultType>(expression.getFirstOperand()->accept(*this, boost::any())).expr()->toExpression() == expression.getManager().integer(index);
                                expressions[index] = std::move(isCurrentIndex);
                            }
                        } else {
                            expressions[expression.getFirstOperand()->evaluateAsInt()] = expression.getManager().boolean(true);
                        }
                    } else {
                        if (expression.getFirstOperand()->containsVariables()) {
                            for (size_t index = 0; index < size; ++index) {
                                storm::expressions::Expression isCurrentIndex = boost::any_cast<ResultType>(expression.getFirstOperand()->accept(*this, boost::any())).expr()->toExpression() == expression.getManager().integer(index);
                                for (auto& entry : resultSecond->first) {
                                    expressions[index * sizeSoFar + entry.first] = (isCurrentIndex && entry.second).simplify();
                                }
                            }
                        } else {
                            auto index = expression.getFirstOperand()->evaluateAsInt();
                            for (auto& entry : resultSecond->first) {
                                expressions[index * sizeSoFar + entry.first] = entry.second;
                            }
                        }
                    }

                    if (arrayNumber == 0) {
                        assert (expressions.size() > 0);
                        assert (expressions.size() > 1 || expressions.begin()->second.evaluateAsBool());
                        auto itr = expressions.begin();
                        storm::expressions::Expression result = boost::any_cast<ResultType>(arrayExpression->accept(*this, itr->first)).expr()->toExpression();
                        ++itr;
                        for (; itr != expressions.end(); ++itr) {
                            result = storm::expressions::ite(itr->second, boost::any_cast<ResultType>(arrayExpression->accept(*this, itr->first)).expr()->toExpression(), result);
                        }
                        return ResultType(result.simplify().getBaseExpressionPointer());
                    } else {
                        std::pair<std::map<uint_fast64_t, storm::expressions::Expression>, uint_fast64_t> result = {expressions, size * sizeSoFar};
                        return result;
                    }
                }
                
                virtual boost::any visit(storm::expressions::FunctionCallExpression const&, boost::any const&) override {
                    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "Found Function call expression while eliminating array expressions. This is not expected since functions are expected to be eliminated at this point.");
                    return false;
                }
                
            private:
                std::unordered_map<storm::expressions::Variable, std::vector<storm::jani::Variable const*>> const& replacements;
                std::unordered_map<storm::expressions::Variable, std::size_t> const& arraySizes;
            };
            
            class MaxArraySizeDeterminer : public ConstJaniTraverser {
            public:
                
                typedef std::unordered_map<storm::expressions::Variable, std::size_t>* MapPtr;
                
                MaxArraySizeDeterminer() = default;
                virtual ~MaxArraySizeDeterminer() = default;
                std::unordered_map<storm::expressions::Variable, std::size_t> getMaxSizes(Model const& model) {
                    // We repeatedly determine the max array sizes until convergence. This is to cover assignments of one array variable to another (A := B)
                    std::unordered_map<storm::expressions::Variable, std::size_t> result, previousResult;
                    do {
                        previousResult = result;
                        ConstJaniTraverser::traverse(model, &result);
                    } while (previousResult != result);
                    return result;
                }
                
                virtual void traverse(Assignment const& assignment, boost::any const& data) override {
                    if (assignment.lValueIsVariable() && assignment.getExpressionVariable().getType().isArrayType()) {
                        auto& map = *boost::any_cast<MapPtr>(data);
                        std::size_t newSize = MaxArraySizeExpressionVisitor().getMaxSize(assignment.getAssignedExpression(), map);
                        auto insertionRes = map.emplace(assignment.getExpressionVariable(), newSize);
                        if (!insertionRes.second) {
                            insertionRes.first->second = std::max(newSize, insertionRes.first->second);
                        }
                    }
                }
                
                virtual void traverse(Variable const& variable, boost::any const& data) override {
                    if (variable.isArrayVariable() && variable.hasInitExpression()) {
                        auto& map = *boost::any_cast<MapPtr>(data);
                        std::size_t newSize = MaxArraySizeExpressionVisitor().getMaxSize(variable.getInitExpression(), map);
                        auto insertionRes = map.emplace(variable.getExpressionVariable(), newSize);
                        if (!insertionRes.second) {
                            insertionRes.first->second = std::max(newSize, insertionRes.first->second);
                        }
                    }
                }
            };
            
        
            class ArrayVariableReplacer : public JaniTraverser {
            public:

                typedef ArrayEliminatorData ResultType;
                using JaniTraverser::traverse;
                
                ArrayVariableReplacer(storm::expressions::ExpressionManager& expressionManager, bool keepNonTrivialArrayAccess, std::unordered_map<storm::expressions::Variable, std::size_t> const& arrayVarToSizesMap) : expressionManager(expressionManager) , keepNonTrivialArrayAccess(keepNonTrivialArrayAccess), arraySizes(arrayVarToSizesMap) {}
                
                virtual ~ArrayVariableReplacer() = default;
                ResultType replace(Model& model) {

                    ResultType result;
                    arrayExprEliminator = std::make_unique<ArrayExpressionEliminationVisitor>(result.replacements, arraySizes);
                    for (auto const& arraySize : arraySizes) {
                        result.replacements.emplace(arraySize.first, std::vector<storm::jani::Variable const*>(arraySize.second, nullptr));
                    }
                    traverse(model, &result);
                    return result;
                }

                virtual void traverse(Model& model, boost::any const& data) override {
                    
                    // Insert fresh basic variables for global array variables
                    auto& replacements = boost::any_cast<ResultType*>(data)->replacements;
                    for (storm::jani::Variable const& arrayVariable : model.getGlobalVariables().getArrayVariables()) {
                        std::vector<storm::jani::Variable const*>& basicVars = replacements.at(arrayVariable.getExpressionVariable());
                        for (uint64_t index = 0; index < basicVars.size(); ++index) {
                            basicVars[index] = &model.addVariable(*getBasicVariable(arrayVariable, index));
                        }
                    }
                    // drop all occuring array variables
                    auto elVars = model.getGlobalVariables().dropAllArrayVariables();
                    auto& eliminatedArrayVariables = boost::any_cast<ResultType*>(data)->eliminatedArrayVariables;
                    eliminatedArrayVariables.insert(eliminatedArrayVariables.end(), elVars.begin(), elVars.end());
                    
                    // Make new variable replacements known to the expression eliminator
                    arrayExprEliminator = std::make_unique<ArrayExpressionEliminationVisitor>(replacements, arraySizes);
                    
                    for (auto& aut : model.getAutomata()) {
                        traverse(aut, data);
                    }
                    
                    // traversal of remaining components
                    if (model.hasInitialStatesRestriction()) {
                        model.setInitialStatesRestriction(arrayExprEliminator->eliminate(model.getInitialStatesRestriction()));
                    }
                    for (auto& nonTrivRew : model.getNonTrivialRewardExpressions()) {
                        nonTrivRew.second = arrayExprEliminator->eliminate(nonTrivRew.second);
                    }
                }
                
                virtual void traverse(Automaton& automaton, boost::any const& data) override {
                    // No need to traverse the init restriction.
                    
                    // Insert fresh basic variables for local array variables
                    auto& replacements = boost::any_cast<ResultType*>(data)->replacements;
                    for (storm::jani::Variable const& arrayVariable : automaton.getVariables().getArrayVariables()) {
                        std::vector<storm::jani::Variable const*>& basicVars = replacements.at(arrayVariable.getExpressionVariable());
                        for (uint64_t index = 0; index < basicVars.size(); ++index) {
                            basicVars[index] = &automaton.addVariable(*getBasicVariable(arrayVariable, index));
                        }
                    }
                    // drop all occuring array variables
                    auto elVars = automaton.getVariables().dropAllArrayVariables();
                    auto& eliminatedArrayVariables = boost::any_cast<ResultType*>(data)->eliminatedArrayVariables;
                    eliminatedArrayVariables.insert(eliminatedArrayVariables.end(), elVars.begin(), elVars.end());

                    // Make new variable replacements known to the expression eliminator
                    arrayExprEliminator = std::make_unique<ArrayExpressionEliminationVisitor>(replacements, arraySizes);
                    
                    for (auto& loc : automaton.getLocations()) {
                        traverse(loc, data);
                    }
                    traverse(automaton.getEdgeContainer(), data);
                    
                    if (automaton.hasInitialStatesRestriction()) {
                        automaton.setInitialStatesRestriction(arrayExprEliminator->eliminate(automaton.getInitialStatesRestriction()));
                    }
                }

                virtual void traverse(Location& location, boost::any const& data) override {
                    traverse(location.getAssignments(), data);
                    if (location.hasTimeProgressInvariant()) {
                        location.setTimeProgressInvariant(arrayExprEliminator->eliminate(location.getTimeProgressInvariant()));
                        traverse(location.getTimeProgressInvariant(), data);
                    }
                }

                void traverse(TemplateEdge& templateEdge, boost::any const& data) override {
                    templateEdge.setGuard(arrayExprEliminator->eliminate(templateEdge.getGuard()));
                    for (auto& dest : templateEdge.getDestinations()) {
                        traverse(dest, data);
                    }
                    traverse(templateEdge.getAssignments(), data);
                }

                
                void traverse(Edge& edge, boost::any const& data) override {
                    if (edge.hasRate()) {
                        edge.setRate(arrayExprEliminator->eliminate(edge.getRate()));
                    }
                    for (auto& dest : edge.getDestinations()) {
                        traverse(dest, data);
                    }
                }
                
                void traverse(EdgeDestination& edgeDestination, boost::any const&) override {
                    edgeDestination.setProbability(arrayExprEliminator->eliminate(edgeDestination.getProbability()));
                }
                
                virtual void traverse(OrderedAssignments& orderedAssignments, boost::any const& data) override {
                    auto const& replacements = boost::any_cast<ResultType*>(data)->replacements;

                    // Replace array occurrences in LValues and assigned expressions.
                    std::vector<Assignment> newAssignments;
                    if (!orderedAssignments.empty()) {
                        int64_t level = orderedAssignments.getLowestLevel();
                        std::unordered_map<storm::expressions::Variable, std::vector<Assignment const*>> collectedArrayAccessAssignments;
                        for (Assignment const& assignment : orderedAssignments) {
                            if (assignment.getLevel() != level) {
                                STORM_LOG_ASSERT(assignment.getLevel() > level, "Ordered Assignment does not have the expected order.");
                                for (auto const& arrayAssignments : collectedArrayAccessAssignments) {
                                    insertLValueArrayAccessReplacements(arrayAssignments.second, replacements.at(arrayAssignments.first), level, newAssignments);
                                }
                                collectedArrayAccessAssignments.clear();
                                level = assignment.getLevel();
                            }
                            if (assignment.getVariable().isArrayVariable() && assignment.getLValue().isFullArrayAccess()) {
                                // This is the easy case, we can eliminate the variable by replacing it by its _at_... equivalent
                                if (!keepNonTrivialArrayAccess || !assignment.getLValue().arrayIndexContainsVariable()) {
                                    auto insertionRes = collectedArrayAccessAssignments.emplace(assignment.getVariable().getExpressionVariable(), std::vector<Assignment const*>({&assignment}));
                                    if (!insertionRes.second) {
                                        insertionRes.first->second.push_back(&assignment);
                                    }
                                } else {
                                    // Keeping array access LValue
                                    auto newIndex = arrayExprEliminator->eliminate(assignment.getLValue().getArrayIndex());
                                    LValue newLValue(LValue(assignment.getVariable(), newIndex, assignment.getLValue().getTotalSize()));
                                    newAssignments.emplace_back(newLValue, arrayExprEliminator->eliminate(assignment.getAssignedExpression()), assignment.getLevel());
                                }
                            } else if (assignment.getVariable().isArrayVariable()) {
                                // In this case we assign to an array to an array, there are one or more indices missing
                                STORM_LOG_ASSERT(assignment.getAssignedExpression().getType().isArrayType(), "Assigning a non-array expression to an array variable... " << assignment);
                                // We can only deal with the case where the last index missing
                                STORM_LOG_THROW((assignment.getLValue().getSizes().size() - 1 == assignment.getLValue().getArrayIndexVector().size()), storm::exceptions::NotImplementedException, "Eliminating nested arrays assigned to nested arrays not yet implemented problem occured at "<< assignment.getName());
                                std::vector<storm::jani::Variable const*> const& arrayVariableReplacements = replacements.at(assignment.getExpressionVariable());
                                // Get the maximum size of the array expression on the rhs
                                uint64_t rhsSize = MaxArraySizeExpressionVisitor().getMaxSize(assignment.getAssignedExpression(), arraySizes);
                                STORM_LOG_ASSERT(arrayVariableReplacements.size() >= rhsSize, "Array size too small.");
                                for (uint64_t index = 0; index < arrayVariableReplacements.size(); ++index) {
                                    auto const& replacement = *arrayVariableReplacements[index];
                                    storm::expressions::Expression newRhs;
                                    if (index < rhsSize) {
                                        newRhs = std::make_shared<storm::expressions::ArrayAccessExpression>(expressionManager, assignment.getAssignedExpression().getType().getElementType(), assignment.getAssignedExpression().getBaseExpressionPointer(), expressionManager.integer(index).getBaseExpressionPointer())->toExpression();
                                    } else {
                                        newRhs = getOutOfBoundsValue(replacement);
                                    }
                                    newRhs = arrayExprEliminator->eliminate(newRhs);
                                    newAssignments.emplace_back(LValue(replacement), newRhs, level);
                                }
                            } else {
                                assert (!assignment.getVariable().isArrayVariable());
                                // In this case we have a non-array variable, we don't need to eliminate anything in the lhs
                                newAssignments.emplace_back(assignment.getLValue(), arrayExprEliminator->eliminate(assignment.getAssignedExpression()), assignment.getLevel());
                            }
                        }
                        for (auto const& arrayAssignments : collectedArrayAccessAssignments) {
                            insertLValueArrayAccessReplacements(arrayAssignments.second, replacements.at(arrayAssignments.first), level, newAssignments);
                        }
                        collectedArrayAccessAssignments.clear();
                        orderedAssignments.clear();
                        for (auto const& assignment : newAssignments) {
                            orderedAssignments.add(assignment);
                            STORM_LOG_ASSERT(!containsArrayExpression(assignment.getExpressionVariable()), "HELP: " << assignment);
                        }
                    }
                }
                
            private:
                
                std::shared_ptr<Variable> getBasicVariable(Variable const& arrayVariable, uint64_t index) const {
                    std::string name = arrayVariable.getExpressionVariable().getName() + "_at_" + std::to_string(index);
                    storm::expressions::Expression initValue;
                    if (arrayVariable.hasInitExpression()) {
                        // We want to eliminate the initial value for the arrayVariable at entry index, in case of a nested_array, this index is already re-calculated.
                        // In the end we have an ArrayAccessExpression
                        auto indexExpression = std::make_shared<storm::expressions::ArrayAccessIndexExpression>(expressionManager, expressionManager.getIntegerType(), expressionManager.integer(index).getBaseExpressionPointer());
                        initValue = arrayExprEliminator->eliminate(std::make_shared<storm::expressions::ArrayAccessExpression>(expressionManager, arrayVariable.getExpressionVariable().getType().getElementType(), arrayVariable.getInitExpression().getBaseExpressionPointer(), indexExpression)->toExpression());
                        STORM_LOG_ASSERT(!containsArrayExpression(initValue), "HELP:  " << arrayVariable.getName());
                    }

                    if (arrayVariable.getType()->getChildType()->isIntegerType()) {
                        storm::expressions::Variable exprVariable = expressionManager.declareIntegerVariable(name);
                        if (arrayVariable.isBoundedVariable()) {
                            return storm::jani::Variable::makeBoundedVariable(name, storm::jani::JaniType::ElementType::Int, exprVariable, initValue, arrayVariable.isTransient(), arrayVariable.getLowerBound(), arrayVariable.getUpperBound());
                        } else {
                            return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Int, exprVariable, initValue, arrayVariable.isTransient());
                        }
                    } else if (arrayVariable.getType()->getChildType()->isRealType()) {
                        storm::expressions::Variable exprVariable = expressionManager.declareRationalVariable(name);
                        return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Real, exprVariable, initValue, arrayVariable.isTransient());
                    } else if (arrayVariable.getType()->getChildType()->isBooleanType()) {
                        storm::expressions::Variable exprVariable = expressionManager.declareBooleanVariable(name);
                        return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Bool, exprVariable, initValue, arrayVariable.isTransient());
                    } else if (arrayVariable.getType()->getChildType()->isArrayType()) {
                        auto const childType = arrayVariable.getType()->getChildType();
                        if (childType->getChildType()->isIntegerType()) {
                            storm::expressions::Variable exprVariable = expressionManager.declareIntegerVariable(name);
                            if (arrayVariable.isBoundedVariable()) {
                                return storm::jani::Variable::makeBoundedVariable(name, storm::jani::JaniType::ElementType::Int, exprVariable, initValue, arrayVariable.isTransient(), arrayVariable.getLowerBound(), arrayVariable.getUpperBound());
                            } else {
                                return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Int, exprVariable, initValue, arrayVariable.isTransient());
                            }
                        } else if (childType->getChildType()->isRealType()) {
                            storm::expressions::Variable exprVariable = expressionManager.declareRationalVariable(name);
                            return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Real, exprVariable, initValue, arrayVariable.isTransient());
                        } else if (childType->getChildType()->isBooleanType()) {
                            storm::expressions::Variable exprVariable = expressionManager.declareBooleanVariable(name);
                            return storm::jani::Variable::makeBasicVariable(name, storm::jani::JaniType::ElementType::Bool, exprVariable, initValue, arrayVariable.isTransient());
                        } else if (childType->getChildType()->isArrayType()) {
                            STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "More than two nested arrays not implemented");
                        }
                    }
                    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Unhandled array base type.");
                    return nullptr;
                }
                
                void insertLValueArrayAccessReplacements(std::vector<Assignment const*> const& arrayAccesses, std::vector<storm::jani::Variable const*> const& arrayVariableReplacements, int64_t level, std::vector<Assignment>& newAssignments) const {
                    bool onlyConstantIndices = true;
                    for (auto const& aa : arrayAccesses) {
                        if (aa->getLValue().arrayIndexContainsVariable()) {
                            onlyConstantIndices = false;
                            break;
                        }
                    }
                    if (onlyConstantIndices) {
                        for (auto const& aa : arrayAccesses) {
                            LValue lvalue(*arrayVariableReplacements.at(aa->getLValue().getArrayIndex().evaluateAsInt()));
                            newAssignments.emplace_back(lvalue, arrayExprEliminator->eliminate(aa->getAssignedExpression()), level);
                        }
                    } else {
                        for (uint64_t index = 0; index < arrayVariableReplacements.size(); ++index) {
                            storm::expressions::Expression assignedExpression = arrayVariableReplacements[index]->getExpressionVariable().getExpression();
                            auto indexExpression = expressionManager.integer(index);
                            for (auto const& aa : arrayAccesses) {
                                assert (false);
                                // TODO: fix this
//                                assignedExpression = storm::expressions::ite(arrayExprEliminator->eliminate(aa->getLValue().getArrayIndex()) == indexExpression, arrayExprEliminator->eliminate(aa->getAssignedExpression()), assignedExpression);
                            }
                            assert (false);
//                            newAssignments.emplace_back(LValue(*arrayVariableReplacements[index]), assignedExpression, level);
                        }
                    }
                }
                
                storm::expressions::Expression getOutOfBoundsValue(Variable const& var) const {
                    if (var.hasInitExpression()) {
                        STORM_LOG_ASSERT(!containsArrayExpression(var.getInitExpression()), "HELP out of bounds " << var.getName());
                        return var.getInitExpression();
                    }
                    if (var.isBooleanVariable()) {
                        return expressionManager.boolean(false);
                    }
                    if (var.isBoundedVariable() && var.isIntegerVariable()) {
                        return var.getLowerBound();
                    }
                    if (!var.isBoundedVariable() && var.isIntegerVariable()) {
                        return expressionManager.integer(0);
                    }
                    if (var.isRealVariable()) {
                        return expressionManager.rational(0.0);
                    }
                    STORM_LOG_THROW(false, storm::exceptions::UnexpectedException, "unhandled variabe type");
                    return storm::expressions::Expression();
                }
                
                std::unique_ptr<ArrayExpressionEliminationVisitor> arrayExprEliminator;
                storm::expressions::ExpressionManager& expressionManager;
                bool const keepNonTrivialArrayAccess;
                std::unordered_map<storm::expressions::Variable, std::size_t> const& arraySizes;
            };
        } // namespace detail
        
        storm::expressions::Expression ArrayEliminatorData::transformExpression(storm::expressions::Expression const& arrayExpression) const {
            std::unordered_map<storm::expressions::Variable, std::size_t> arraySizes;
            for (auto const& r : replacements) {
                arraySizes.emplace(r.first, r.second.size());
            }
            detail::ArrayExpressionEliminationVisitor eliminator(replacements, arraySizes);
            return eliminator.eliminate(arrayExpression);
        }
        
        void ArrayEliminatorData::transformProperty(storm::jani::Property& property) const {
            property = property.substitute([this](storm::expressions::Expression const& exp) {return transformExpression(exp);});
        }
        
        ArrayEliminatorData ArrayEliminator::eliminate(Model& model, bool keepNonTrivialArrayAccess) {
            ArrayEliminatorData result;
            // Only perform actions if there actually are arrays.
            if (model.getModelFeatures().hasArrays()) {
                auto sizes = detail::MaxArraySizeDeterminer().getMaxSizes(model);
                result = detail::ArrayVariableReplacer(model.getExpressionManager(), keepNonTrivialArrayAccess, sizes).replace(model);
                if (!keepNonTrivialArrayAccess) {
                    model.getModelFeatures().remove(ModelFeature::Arrays);
                }
                // TODO: hack
                auto size = model.getConstants().size();
                for (size_t i = 0; i < size; ++i) {
                    auto constant = model.getConstants().at(i);
                    if (constant.isDefined() && containsArrayExpression(constant.getExpression())) {
                        // We hope we don't need this one any longer however this could break everything :D Breaking things is fun :D
                        model.removeConstant(constant.getName());
                        i--;
                        size--;
                    }
                }
                model.finalize();
            }
            STORM_LOG_ASSERT(!containsArrayExpression(model), "the model still contains array expressions.");
            return result;
        }
    }
}

