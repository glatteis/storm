#ifndef STORM_STORAGE_PRISM_BOOLEANVARIABLE_H_
#define STORM_STORAGE_PRISM_BOOLEANVARIABLE_H_

#include <map>

#include "src/storage/prism/Variable.h"

namespace storm {
    namespace prism {
        class BooleanVariable : public Variable {
        public:
            // Create default implementations of constructors/assignment.
            BooleanVariable() = default;
            BooleanVariable(BooleanVariable const& other) = default;
            BooleanVariable& operator=(BooleanVariable const& other)= default;
            BooleanVariable(BooleanVariable&& other) = default;
            BooleanVariable& operator=(BooleanVariable&& other) = default;

            /*!
             * Creates a boolean variable with the given name and the default initial value expression.
             *
             * @param variableName The name of the variable.
             * @param filename The filename in which the variable is defined.
             * @param lineNumber The line number in which the variable is defined.
             */
            BooleanVariable(std::string const& variableName, std::string const& filename = "", uint_fast64_t lineNumber = 0);

            /*!
             * Creates a boolean variable with the given name and the given constant initial value expression.
             *
             * @param variableName The name of the variable.
             * @param initialValueExpression The constant expression that defines the initial value of the variable.
             * @param filename The filename in which the variable is defined.
             * @param lineNumber The line number in which the variable is defined.
             */
            BooleanVariable(std::string const& variableName, storm::expressions::Expression const& initialValueExpression, std::string const& filename = "", uint_fast64_t lineNumber = 0);
            
            /*!
             * Creates a copy of the given boolean variable and performs the provided renaming.
             *
             * @param oldVariable The variable to copy.
             * @param newName New name of this variable.
             * @param renaming A mapping from names that are to be renamed to the names they are to be
             * replaced with.
             * @param filename The filename in which the variable is defined.
             * @param lineNumber The line number in which the variable is defined.
             */
            BooleanVariable(BooleanVariable const& oldVariable, std::string const& newName, std::map<std::string, std::string> const& renaming, std::string const& filename, uint_fast64_t lineNumber);
            
            friend std::ostream& operator<<(std::ostream& stream, BooleanVariable const& variable);
        };
        
    } // namespace prism
} // namespace storm

#endif /* STORM_STORAGE_PRISM_BOOLEANVARIABLE_H_ */
