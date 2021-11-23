#include "test/storm_gtest.h"
#include "logic/Formula.h"
#include "models/sparse/Dtmc.h"
#include "storage/prism/Program.h"
#include "storm-config.h"
#include "test/storm_gtest.h"

#include "storm-pars/api/analysis.h"
#include "storm-pars/api/region.h"
#include "storm-pars/transformer/SparseParametricDtmcSimplifier.h"

#include "storm-parsers/api/storm-parsers.h"
#include "storm-parsers/parser/AutoParser.h"
#include "storm-parsers/parser/PrismParser.h"

#include "utility/prism.h"
#include "storm-pars/transformer/RationalTransitionUnfolder.h"
namespace {

		TEST(RationalTransitionUnfolderTest, TestRationalFunctionUnfolding) {
				std::string programFile = STORM_TEST_RESOURCES_DIR "/pdtmc/brp16_2.pm";
				std::string formulaAsString = "P=? [F s=4 & i=N ]";
				std::string constantsAsString = ""; //e.g. pL=0.9,TOACK=0.5

				// Program and formula
				storm::prism::Program program = storm::api::parseProgram(programFile);
				program = storm::utility::prism::preprocess(program, constantsAsString);
				std::vector<std::shared_ptr<const storm::logic::Formula>> formulas = storm::api::extractFormulasFromProperties(storm::api::parsePropertiesForPrismProgram(formulaAsString, program));
				std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> model = storm::api::buildSparseModel<storm::RationalFunction>(program, formulas)->as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
				std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> dtmc = model->as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
				auto simplifier = storm::transformer::SparseParametricDtmcSimplifier<storm::models::sparse::Dtmc<storm::RationalFunction>>(*dtmc);
				ASSERT_TRUE(simplifier.simplify(*(formulas[0])));
				model = simplifier.getSimplifiedModel();

				dtmc = model->as<storm::models::sparse::Dtmc<storm::RationalFunction>>();

				ASSERT_EQ(193ul, dtmc->getNumberOfStates());
				ASSERT_EQ(383ul, dtmc->getNumberOfTransitions());
				storm::transformer::RationalTransitionUnfolder<storm::models::sparse::Dtmc<storm::RationalFunction>> transitionUnfolder(*dtmc);

				std::shared_ptr<storm::RawPolynomialCache> cache = std::make_shared<storm::RawPolynomialCache>();
				auto p = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(carl::VariablePool::getInstance().getFreshPersistentVariable("vp")), cache));
				auto q = storm::RationalFunction(storm::Polynomial(storm::RawPolynomial(carl::VariablePool::getInstance().getFreshPersistentVariable("vq")), cache));
				/* auto polynomial = storm::utility::convertNumber<storm::RationalFunction>(p * p * (storm::utility::one<storm::RationalFunction>() - q) * 3); */
				storm::parser::ValueParser<storm::RationalFunction> functionParser;
				functionParser.addParameter("p");
				storm::RationalFunction polynomial = functionParser.parseValue("(5 * ((p)^2 * (p+(-1))^2))/(1)");

				std::cout << "unfolding: " << polynomial << std::endl;;
				auto result = transitionUnfolder.unfoldIntoSimpleChunks(polynomial);
				std::cout << "is initialized: " << result.is_initialized() << std::endl;

				for (auto const& member : *result) {
					std::cout << "member: " << member << std::endl;
				}
		}
}
