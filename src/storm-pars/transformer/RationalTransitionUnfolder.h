#pragma once
#include <memory>
#include "adapters/RationalFunctionAdapter.h"
namespace storm {
    namespace transformer {

        /*!
         * The RationalTransitionUnfolder unfolds transitions that have rational functions into states with constant transitions, transitions of the form p or (1-p).
         * This fails if the rational functions cannot be broken apart.
         */
        template<typename SparseModelType>
        class RationalTransitionUnfolder {
        public:
            RationalTransitionUnfolder(SparseModelType const& model);
            boost::optional<std::vector<storm::RationalFunction>> unfoldIntoSimpleChunks(storm::RationalFunction const& f);
        private:
            SparseModelType const& originalModel;
            std::shared_ptr<SparseModelType> unfoldedModel;
        };
		}
}
