#include "storm-pars/transformer/RationalTransitionUnfolder.h"
#include "models/sparse/Dtmc.h"
#include "models/sparse/Mdp.h"
#include "utility/constants.h"
#include "utility/logging.h"

namespace storm {
    namespace transformer {
        template<typename SparseModelType>
        RationalTransitionUnfolder<SparseModelType>::RationalTransitionUnfolder(SparseModelType const& model) : originalModel(model) {
            // intentionally left empty

        }

        template<typename SparseModelType>
        boost::optional<std::vector<storm::RationalFunction>> RationalTransitionUnfolder<SparseModelType>::unfoldIntoSimpleChunks(storm::RationalFunction const& f) {
            if (!f.denominatorAsPolynomial().isConstant()) {
                STORM_LOG_INFO("Could not unfold transitions as rational function is not a polynomial.");
                return boost::none;
            }
            std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
            RawPolynomial workingFormula(f.nominatorAsPolynomial());
            std::vector<storm::RationalFunction> extractedChunks;
            while (true) {
                if (workingFormula.isConstant()) {
                    break;
                }
                auto variables = workingFormula.gatherVariables();
                auto nextVariable = storm::RawPolynomial(*variables.begin());

                bool successInUnfolding = false;

                std::vector<storm::RawPolynomial> functionsToDivideBy = {nextVariable, RawPolynomial(1) - nextVariable};
                for (auto const& functionToDivideBy : functionsToDivideBy) {
                    carl::DivisionResult<RawPolynomial> divisionResultVariable = workingFormula.divideBy(functionToDivideBy);
                    if (divisionResultVariable.remainder.isZero()) {
                        successInUnfolding = true;
                        extractedChunks.push_back(storm::RationalFunction(storm::Polynomial(functionToDivideBy, cache)));
                        workingFormula = divisionResultVariable.quotient;
                    }
                }

                if (!successInUnfolding) {
                    STORM_LOG_INFO("Could not unfold variable " << nextVariable);
                    return boost::none;
                }
            }
            // Add constants if nessecary
            if (!workingFormula.isOne()) {
                extractedChunks.push_back(storm::RationalFunction(storm::Polynomial(workingFormula, cache)));
            }
            if (!f.denominatorAsPolynomial().isOne()) {
                extractedChunks.push_back(storm::RationalFunction(f.denominatorAsPolynomial()));
            }
            return extractedChunks;
        }
        template class RationalTransitionUnfolder<models::sparse::Dtmc<storm::RationalFunction>>;
        template class RationalTransitionUnfolder<models::sparse::Mdp<storm::RationalFunction>>;
		}
}
