#ifndef STORM_EQUALPARAMATERREDUCER_H
#define STORM_EQUALPARAMATERREDUCER_H

#include <set>
#include "adapters/RationalFunctionAdapter.h"
#include "modelchecker/CheckTask.h"
#include "models/sparse/Dtmc.h"
#include "storage/FlexibleSparseMatrix.h"
#include "storm-pars/utility/parametric.h"

namespace storm {
namespace transformer {

class TimeTravelling {
    
   public:
        TimeTravelling() {}
        models::sparse::Dtmc<RationalFunction> timeTravel(models::sparse::Dtmc<RationalFunction> model,
                                                          modelchecker::CheckTask<logic::Formula, RationalNumber> const& checkTask);

   private:
        void updateTreeStates(std::map<RationalFunctionVariable, std::map<uint_fast64_t, std::set<uint_fast64_t>>>& treeStates,
                              std::map<RationalFunctionVariable, std::set<uint_fast64_t>>& workingSets,
                              storage::FlexibleSparseMatrix<RationalFunction>& flexibleMatrix, const std::set<carl::Variable>& allParameters,
                              const boost::optional<std::vector<RationalFunction>>& stateRewardVector, const models::sparse::StateLabeling stateLabelling,
                              const std::set<std::string> labelsInFormula);
        std::vector<storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> joinDuplicateTransitions(std::vector<storm::storage::MatrixEntry<uint_fast64_t, RationalFunction>> const& entries);
        bool jipConvert(uint_fast64_t state, storage::FlexibleSparseMatrix<RationalFunction>& matrix, std::map<uint_fast64_t, bool>& alreadyVisited,
                        const std::map<RationalFunctionVariable, std::map<uint_fast64_t, std::set<uint_fast64_t>>>& treeStates,
                        const std::set<carl::Variable>& allParameters, const boost::optional<std::vector<RationalFunction>>& stateRewardVector,
                        const models::sparse::StateLabeling stateLabelling, const std::set<std::string> labelsInFormula);
        bool labelsIntersectedEqual(const std::set<std::string>& labels1, const std::set<std::string>& labels2, const std::set<std::string>& intersection);
};

}  // namespace transformer
}  // namespace storm
    
#endif
