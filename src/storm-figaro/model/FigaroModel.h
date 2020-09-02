#pragma once
<<<<<<< HEAD
#include "FigaroModelTemplate.h"
=======

>>>>>>> 2dddc29ea061c2a0b33493dccd8ece2b3a1211bf
#include <array>
#include <map>
#include <vector>
#include <sstream>
#include <math.h>
#include <set>

<<<<<<< HEAD
namespace storm{
    namespace figaro{
        class FigaroProgram1: public storm::figaro::FigaroProgram{
        public:
            FigaroProgram1(): FigaroProgram(


//            std::map<std::string, size_t> mFigaroboolelementindex =
                    {},

//            std::map<std::string, size_t> mFigaroelementfailureindex =
                    { {"mark_OF_Second_waiting_area_iseq_3",0}},

//            std::map<std::string, size_t> mFigarofloatelementindex =
                    {
                            {"calculated_lambda_OF_Arrival" , 0},
                            {"calculated_lambda_OF_end_of_payment" , 1},
                            {"calculated_lambda_OF_tank_is_full_1" , 2},
                            {"calculated_lambda_OF_tank_is_full_2" , 3},
                            {"calculated_lambda_OF_tank_is_full_3" , 4}},

//            std::map<std::string, size_t> mFigarointelementindex =
                    {
                            {"mark_OF_Cashdesk" , 0},
                            {"mark_OF_First_waiting_area" , 1},
                            {"mark_OF_Pump_1" , 2},
                            {"mark_OF_Pump_2" , 3},
                            {"mark_OF_Pump_3" , 4},
                            {"mark_OF_Second_waiting_area" , 5}},

//            std::map<std::string, size_t> mFigaroenumelementindex =
                    { },
//            std::map<std::string, size_t> failure_variable_names =

                    {"mark_OF_Second_waiting_area_iseq_3"},
//            std::set<std::string> enum_variables_names =
                    { },

//            std::set<std::string> float_variables_names =
                    {
                            "calculated_lambda_OF_Arrival" ,
                            "calculated_lambda_OF_end_of_payment" ,
                            "calculated_lambda_OF_tank_is_full_1" ,
                            "calculated_lambda_OF_tank_is_full_2" ,
                            "calculated_lambda_OF_tank_is_full_3" },
//                std::string const topevent=
                    "mark_OF_Second_waiting_area_iseq_3",
//                static int const
//                        numBoolState =
                            0 ,
//                    numBoolFailureState =
                            1,
//                static int const
//                    numFloatState =
                    5 ,
//                static int const
//                    numIntState =
                    6 ,
//                static int const
//                    numEnumState =
                    0 ,

//            bool ins_transition_found =
                    false){}

            /* ---------- CODING ENUMERATED VARIABLES STATES ------------ */
            enum enum_status {};

//
//            std::array<bool, numBoolState> boolState;
//            std::array<bool, numBoolState> backupBoolState;
//            std::array<double, numFloatState> floatState;
//            std::array<double, numFloatState> backupFloatState;
//            std::array<int, numIntState> intState;
//            std::array<int, numIntState> backupIntState;
//            std::array<int, numEnumState> enumState;
//           ;


=======
#include "FigaroModelTemplate.h"

namespace storm {
    namespace figaro {
        
        class FigaroProgram1: public storm::figaro::FigaroProgram{
        public:
            // Constructor
            FigaroProgram1(): FigaroProgram(
                {
                    {"fail_OF_Node_1" , 0},
                    {"connected_OF_Node_1" , 1},
                    {"fail_OF_Node_2" , 2},
                    {"connected_OF_Node_2" , 3},
                    {"interruption_OF_ud_1" , 4},
                    {"interruption_OF_bidir_3" , 5},
                    {"fail_OF_Node_6" , 6},
                    {"connected_OF_Node_6" , 7},
                    {"fail_OF_Node_8" , 8},
                    {"connected_OF_Node_8" , 9},
                    {"fail_OF_Source" , 10},
                    {"connected_OF_Source" , 11},
                    {"fail_OF_Target" , 12},
                    {"connected_OF_Target" , 13},
                    {"interruption_OF_ud_3" , 14},
                    {"interruption_OF_bidir_7" , 15},
                    {"interruption_OF_bidir_12" , 16},
                    {"interruption_OF_bidir_2" , 17},
                    {"not_connected_OF_Target", 18}
                },
                {{"not_connected_OF_Target", 18}},
                {},
                {{"nb_failures_OF_Failure_counter" , 0}},
                {},
                {},
                {},
                "not_connected_OF_Target",
                19,
                0,
                1,
                0,
                false
            ) {
                // body left empty
            }


            /* ---------- CODING ENUMERATED VARIABLES STATES ------------ */
            enum enum_status {};
            bool REINITIALISATION_OF_connected_OF_Node_1 ;
            bool REINITIALISATION_OF_connected_OF_Node_2 ;
            bool REINITIALISATION_OF_connected_OF_Node_6 ;
            bool REINITIALISATION_OF_connected_OF_Node_8 ;
            bool REINITIALISATION_OF_connected_OF_Source ;
            bool REINITIALISATION_OF_connected_OF_Target ;
            

>>>>>>> 2dddc29ea061c2a0b33493dccd8ece2b3a1211bf
            /* ---------- DECLARATION OF CONSTANTS ------------ */
            int const weight_OF_ai_6 = 3;
            int const weight_OF_ia_1_1 = 1;
            int const weight_OF_am_6 = 1;
            int const weight_OF_av_4 = 1;
            int const weight_OF_am_7 = 1;
            int const weight_OF_ia_1_2 = 1;
            int const weight_OF_am_8 = 1;
            int const weight_OF_am_2 = 1;
            int const weight_OF_av_1 = 1;
            int const weight_OF_ia_4 = 3;
            double const lambda_OF_tank_is_full_3 = 0.1;
            double const lambda_OF_end_of_payment = 0.1;
            int const weight_OF_ai_6_C2 = 3;
            double const lambda_OF_tank_is_full_1 = 0.1;
            int const weight_OF_av_8 = 1;
            int const weight_OF_ia_5 = 2;
            double const lambda_OF_tank_is_full_2 = 0.1;
            int const weight_OF_am_4 = 1;
            int const weight_OF_av_2 = 1;
            int const weight_OF_av_5 = 1;
            int const weight_OF_av_6 = 1;
            int const weight_OF_am_5 = 1;
            bool const Profil1_OF___ARBRE__EIRM = true;
<<<<<<< HEAD
            int const weight_OF_am_3 = 1;
            int const weight_OF_am_1 = 1;
            int const weight_OF_ia_1 = 1;
            int const weight_OF_av_3 = 1;
            int const weight_OF_av_7 = 1;
            int const weight_OF_ai_6_C1 = 3;
            double const lambda_OF_Arrival = 0.2;

            /* ---------- DECLARATION OF OCCURRENCE RULES FIRING FLAGS ------------ */
            bool FIRE_xx2_OF_Arrival;
            bool FIRE_xx2_OF_end_of_payment;
            bool FIRE_xx2_OF_tank_is_full_1;
            bool FIRE_xx2_OF_tank_is_full_2;
            bool FIRE_xx2_OF_tank_is_full_3;


            int calculated_lambda_OF_Arrival = 0 ;
            int calculated_lambda_OF_end_of_payment = 1 ;
            int calculated_lambda_OF_tank_is_full_1 = 2 ;
            int calculated_lambda_OF_tank_is_full_2 = 3 ;
            int calculated_lambda_OF_tank_is_full_3 = 4 ;

            int mark_OF_Cashdesk = 0 ;
            int mark_OF_First_waiting_area = 1 ;
            int mark_OF_Pump_1 = 2 ;
            int mark_OF_Pump_2 = 3 ;
            int mark_OF_Pump_3 = 4 ;
            int mark_OF_Second_waiting_area = 5 ;
            int mark_OF_Second_waiting_area_iseq_3 = 0;


=======
            std::string const function_OF_Node_2 = "intermediate";
            double const link_lambda_OF_ud_1 = 1e-05;
            double const link_mu_OF_bidir_3 = 1;
            double const link_mu_OF_ud_1 = 1;
            double const mu_OF_Node_2 = 0.1;
            double const mu_OF_Target = 0.1;
            double const lambda_OF_Node_1 = 1e-05;
            double const lambda_OF_Source = 1e-05;
            double const lambda_OF_Node_6 = 1e-05;
            double const link_lambda_OF_bidir_2 = 1e-05;
            std::string const function_OF_Target = "target";
            double const link_lambda_OF_ud_3 = 1e-05;
            std::string const function_OF_Node_1 = "intermediate";
            double const lambda_OF_Node_8 = 1e-05;
            std::string const function_OF_Node_6 = "intermediate";
            double const link_lambda_OF_bidir_12 = 1e-05;
            double const link_mu_OF_bidir_2 = 1;
            double const mu_OF_Node_1 = 0.1;
            std::string const function_OF_Source = "source";
            

            /* ---------- DECLARATION OF OCCURRENCE RULES FIRING FLAGS ------------ */
            bool FIRE_xx1_OF_Node_1;
            bool FIRE_xx2_OF_Node_1;
            bool FIRE_xx1_OF_Node_2;
            bool FIRE_xx2_OF_Node_2;
            bool FIRE_xx3_OF_ud_1;
            bool FIRE_xx4_OF_ud_1;
            bool FIRE_xx3_OF_bidir_3;
            bool FIRE_xx4_OF_bidir_3;
            bool FIRE_xx1_OF_Node_6;
            bool FIRE_xx2_OF_Node_6;
            bool FIRE_xx1_OF_Node_8;
            bool FIRE_xx2_OF_Node_8;
            bool FIRE_xx1_OF_Source;
            bool FIRE_xx2_OF_Source;
            bool FIRE_xx1_OF_Target;
            bool FIRE_xx2_OF_Target;
            bool FIRE_xx3_OF_ud_3;
            bool FIRE_xx4_OF_ud_3;
            bool FIRE_xx3_OF_bidir_7;
            bool FIRE_xx4_OF_bidir_7;
            bool FIRE_xx3_OF_bidir_12;
            bool FIRE_xx4_OF_bidir_12;
            bool FIRE_xx3_OF_bidir_2;
            bool FIRE_xx4_OF_bidir_2;
            
            int fail_OF_Node_1 = 0 ;
            int connected_OF_Node_1 = 1 ;
            int fail_OF_Node_2 = 2 ;
            int connected_OF_Node_2 = 3 ;
            int interruption_OF_ud_1 = 4 ;
            int interruption_OF_bidir_3 = 5 ;
            int fail_OF_Node_6 = 6 ;
            int connected_OF_Node_6 = 7 ;
            int fail_OF_Node_8 = 8 ;
            int connected_OF_Node_8 = 9 ;
            int fail_OF_Source = 10 ;
            int connected_OF_Source = 11 ;
            int fail_OF_Target = 12 ;
            int connected_OF_Target = 13 ;
            int interruption_OF_ud_3 = 14 ;
            int interruption_OF_bidir_7 = 15 ;
            int interruption_OF_bidir_12 = 16 ;
            int interruption_OF_bidir_2 = 17 ;
            int not_connected_OF_Target = 18;
            
            int nb_failures_OF_Failure_counter = 0 ;
            
            
>>>>>>> 2dddc29ea061c2a0b33493dccd8ece2b3a1211bf
            /* ---------- DECLARATION OF FUNCTIONS ------------ */
            // Override virtual functions of base class
            void init();
            void saveCurrentState();
            void printState();
            void fireOccurrence(int numFire);
            std::vector<std::tuple<int, double, std::string, int>> showFireableOccurrences();
            void runOnceInteractionStep_default_step();
            int compareStates();
            void doReinitialisations();
            void runInteractions();
            void printstatetuple();
            void fireinsttransitiongroup(std::string);
            int_fast64_t stateSize() const;
            bool figaromodelhasinstransitions();
        };
    }
}

