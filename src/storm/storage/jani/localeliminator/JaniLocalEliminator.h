#pragma once

#include <queue>
#include "storm/storage/jani/Model.h"
#include "storm/storage/jani/Property.h"
#include "boost/variant.hpp"

namespace storm {
    namespace jani {
        class JaniLocalEliminator{
        private:
            class AutomatonInfo{
            public:
                explicit AutomatonInfo();
                std::set<uint64_t> potentiallyPartOfProp;
                bool hasSink;
                uint64_t sinkIndex;
            };
        public:
            class Session {
            public:
                explicit Session(Model model, Property property, bool flatten = true);
                Model &getModel();
                void setModel(const Model &model);
                Property &getProperty();
                bool getFinished() const;
                void setFinished(bool finished);

                bool isLogEnabled();
                void addToLog(const std::string& item);
                std::vector<std::string> getLog();

                AutomatonInfo &getAutomatonInfo(const std::string& name);
                void buildAutomataInfo();
                void flatten_automata();
                void addMissingGuards(const std::string& automatonName);

                expressions::Expression getNewGuard(const Edge& edge, const EdgeDestination& dest, const Edge& outgoing);
                expressions::Expression getProbability(const EdgeDestination& first, const EdgeDestination& then);
                OrderedAssignments executeInSequence(const EdgeDestination& first, const EdgeDestination& then, std::set<std::string> &rewardVariables);
                bool isEliminable(const std::string &automatonName, std::string const& locationName);
                bool hasLoops(const std::string &automatonName, std::string const& locationName);
                bool hasNamedActions(const std::string &automatonName, std::string const& locationName);
                bool isPossiblyInitial(const std::string &automatonName, std::string const &locationName);
                bool isPartOfProp(const std::string &automatonName, std::string const &locationName);
                bool isPartOfProp(const std::string &automatonName, uint64_t locationIndex);
                bool computeIsPartOfProp(const std::string &automatonName, const std::string &locationName);
                bool computeIsPartOfProp(const std::string &automatonName, uint64_t locationIndex);
                bool computeIsPartOfProp(const std::map<expressions::Variable, expressions::Expression>& substitutionMap);
                void setPartOfProp(const std::string &automatonName, const std::string &locationName, bool isPartOfProp);
                void setPartOfProp(const std::string &automatonName, uint64_t locationIndex, bool isPartOfProp);
                void clearIsPartOfProp(const std::string &automatonName);
                bool isVariablePartOfProperty(const std::string &expressionVariableName);

                bool isRewardFormula;
                std::set<std::string> rewardModels;
            private:
                Model model;
                Property property;
                bool finished;
                // We keep a log separate from the main log to prevent the main log from being overwhelmed. This log
                // is exposed via the python API
                std::vector<std::string> log;
                bool logEnabled;

                std::map<std::string, AutomatonInfo> automataInfo;
                std::set<uint_fast64_t> expressionVarsInProperty;
            };

        public:
            class Action {
            public:
                virtual std::string getDescription() = 0;
                virtual void doAction(Session &session) = 0;
            };

            class EliminationScheduler {
            public:
                EliminationScheduler();
                std::unique_ptr<Action> getNextAction();
                void addAction(std::unique_ptr<Action> action);
            private:
                std::queue<std::unique_ptr<Action>> actionQueue;
            };

            EliminationScheduler scheduler;
            explicit JaniLocalEliminator(Model const& original, storm::jani::Property &property, bool addMissingGuards = false);
            explicit JaniLocalEliminator(Model const& original, std::vector<storm::jani::Property>& properties, bool addMissingGuards = false);
            void eliminate(bool flatten = true, bool useTransientVariables = true);
            Model const& getResult();
            std::vector<std::string> getLog();

        private:
            Model const& original;
            Model newModel;
            Property property;
            bool addMissingGuards;
            std::vector<std::string> log;

            void setProperty(storm::jani::Property &property);
        };


    }
}