// EnergyPlus, Copyright (c) 1996-2019, The Board of Trustees of the University of Illinois,
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy), Oak Ridge
// National Laboratory, managed by UT-Battelle, Alliance for Sustainable Energy, LLC, and other
// contributors. All rights reserved.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without the U.S. Department of Energy's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>

// EnergyPlus Headers
#include <DataGlobals.hh>
#include <DataHeatBalance.hh>
#include <DataPrecisionGlobals.hh>
#include <DataRoomAirModel.hh>
#include <General.hh>
#include <HeatBalanceManager.hh>
#include <HybridModel.hh>
#include <InputProcessing/InputProcessor.hh>
#include <OutputProcessor.hh>
#include <ScheduleManager.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

namespace HybridModel {

    // MODULE INFORMATION:
    //       AUTHOR         Sang Hoon Lee, Tianzhen Hong, Rongpeng Zhang. LBNL
    //       DATE WRITTEN   Oct 2015

    // PURPOSE OF THIS MODULE:
    // This module manages hybrid model.

    // METHODOLOGY EMPLOYED:
    //  The model uses measured zone air temperature to calculate internal thermal mass or infiltration air flow rate.

    // USE STATEMENTS:

    // Using/Aliasing
    using namespace DataGlobals;
    using namespace DataHeatBalance;
    using namespace DataPrecisionGlobals;
    using namespace DataRoomAirModel;
    using DataGlobals::ScheduleAlwaysOn;
    using General::CheckCreatedZoneItemName;

    bool FlagHybridModel(false);    // True if hybrid model is activated
    bool FlagHybridModel_TM(false); // User input IM option - True if hybrid model (thermal mass) is activated
    bool FlagHybridModel_AI(false); // User input IM option - True if hybrid model (air infiltration) is activated
    bool FlagHybridModel_PC(false); // User input IM option - True if hybrid model (people count) is activated

    int NumOfHybridModelZones(0);    // Number of hybrid model zones in the model
    std::string CurrentModuleObject; // to assist in getting input

    // SUBROUTINE LOCAL VARIABLE DECLARATIONS:

    // Object Data
    Array1D<HybridModelProperties> HybridModelZone;

    // Functions

    void GetHybridModelZone()
    {

        using ScheduleManager::GetScheduleIndex;

        bool ErrorsFound(false); // If errors detected in input
        Array1D_bool lAlphaFieldBlanks(16, false);
        Array1D_bool lNumericFieldBlanks(4, false);
        int NumAlphas;  // Number of Alphas for each GetobjectItem call
        int NumNumbers; // Number of Numbers for each GetobjectItem call
        int IOStatus;
        int ZonePtr;                     // Pointer to the zone
        int ZoneListPtr;                 // Pointer to the zone list
        std::string CurrentModuleObject; // to assist in getting input
        Array1D_string cAlphaArgs(16);   // Alpha input items for object
        Array1D_string cAlphaFieldNames(16);
        Array1D_string cNumericFieldNames(16);
        Array1D<Real64> rNumericArgs(4); // Numeric input items for object
        int HybridModelStartMonth(0);    // Hybrid model start month
        int HybridModelStartDate(0);     // Hybrid model start date of month
        int HybridModelEndMonth(0);      // Hybrid model end month
        int HybridModelEndDate(0);       // Hybrid model end date of month
        int HMStartDay(0);
        int HMEndDay(0);

        int TemperatureSchPtr(0);      // Temperature schedule pointer
        int HumidityRatioSchPtr(0);    // Humidity ratio schedule pointer
        int CO2ConcentrationSchPtr(0); // CO2 concentration schedule pointer

        int PeopleActivityLevelSchPtr(0);    // People activity level schedule pointer
        int PeopleSensibleFractionSchPtr(0); // People sensible heat portion schedule pointer
        int PeopleRadiantFractionSchPtr(0);  // People radiant heat portion (of sensible heat) schedule pointer
        int PeopleCO2GenRateSchPtr(0);       // People CO2 generation rate schedule pointer

        int SupplyAirTemperatureSchPtr(0);
        int SupplyAirMassFlowRateSchPtr(0);
        int SupplyAirHumidityRatioSchPtr(0);
        int SupplyAirCO2ConcentrationSchPtr(0);

        bool helper_InternalThermalMassCalc(false); // Calculate zone thermal mass flag
        bool helper_AirInfiltrationCalc(false);     // Calculate zone air infiltration rate flag
        bool helper_PeopleCountCalc(false);         // Calculate zone people count flag

        bool helper_IncludeSystemSupply(false); // Include system supply parameters flag

        // Read hybrid model input
        CurrentModuleObject = "HybridModel:Zone";
        NumOfHybridModelZones = inputProcessor->getNumObjectsFound(CurrentModuleObject);
        HybridModelZone.allocate(NumOfZones);

        if (NumOfHybridModelZones > 0) {

            for (int HybridModelNum = 1; HybridModelNum <= NumOfHybridModelZones; ++HybridModelNum) {

                inputProcessor->getObjectItem(CurrentModuleObject,
                                              HybridModelNum,
                                              cAlphaArgs,
                                              NumAlphas,
                                              rNumericArgs,
                                              NumNumbers,
                                              IOStatus,
                                              lNumericFieldBlanks,
                                              lAlphaFieldBlanks,
                                              cAlphaFieldNames,
                                              cNumericFieldNames);

                ZoneListPtr = 0;
                ZonePtr = UtilityRoutines::FindItemInList(cAlphaArgs(2), Zone); // "Zone" is a 1D array, cAlphaArgs(2) is the zone name
                if (ZonePtr == 0 && NumOfZoneLists > 0) ZoneListPtr = UtilityRoutines::FindItemInList(cAlphaArgs(2), ZoneList);
                if (ZonePtr > 0) {
                    HybridModelZone(ZonePtr).Name = cAlphaArgs(1);                          // Zone HybridModel name
                    FlagHybridModel_TM = UtilityRoutines::SameString(cAlphaArgs(3), "Yes"); // Calculate thermal mass option
                    FlagHybridModel_AI = UtilityRoutines::SameString(cAlphaArgs(4), "Yes"); // Calculate infiltration rate option
                    FlagHybridModel_PC = UtilityRoutines::SameString(cAlphaArgs(5), "Yes"); // Calculate people count option

                    // Pointers used to help decide which unknown parameter to solve
                    // Zone Air Infiltration Rate and Zone Internal Thermal Mass calculations cannot be performed simultaneously
                    TemperatureSchPtr = GetScheduleIndex(cAlphaArgs(6));
                    HumidityRatioSchPtr = GetScheduleIndex(cAlphaArgs(7));
                    CO2ConcentrationSchPtr = GetScheduleIndex(cAlphaArgs(8));

                    // Not used for now
                    PeopleActivityLevelSchPtr = GetScheduleIndex(cAlphaArgs(9));
                    PeopleSensibleFractionSchPtr = GetScheduleIndex(cAlphaArgs(10));
                    PeopleRadiantFractionSchPtr = GetScheduleIndex(cAlphaArgs(11));
                    PeopleCO2GenRateSchPtr = GetScheduleIndex(cAlphaArgs(12));

                    // Pointers used to help decide wheather to include system supply terms in the inverse algorithms
                    SupplyAirTemperatureSchPtr = GetScheduleIndex(cAlphaArgs(13));
                    SupplyAirMassFlowRateSchPtr = GetScheduleIndex(cAlphaArgs(14));
                    SupplyAirHumidityRatioSchPtr = GetScheduleIndex(cAlphaArgs(15));
                    SupplyAirCO2ConcentrationSchPtr = GetScheduleIndex(cAlphaArgs(16));

                    /*	Note: Internal thermal mass can be calculated only with measured temperature.
                                      Air infiltration rate can be calculated with either measured temperature, humifity ratio, or CO2 concentration.
                                      People count can be calculated with either measured temperature, humifity ratio, or CO2 concentration.
                    */
                    // Initially set all flags to be false
                    HybridModelZone(ZonePtr).InternalThermalMassCalc_T = false;
                    HybridModelZone(ZonePtr).InfiltrationCalc_T = false;
                    HybridModelZone(ZonePtr).InfiltrationCalc_H = false;
                    HybridModelZone(ZonePtr).InfiltrationCalc_C = false;
                    HybridModelZone(ZonePtr).PeopelCountCalc_T = false;
                    HybridModelZone(ZonePtr).PeopelCountCalc_H = false;
                    HybridModelZone(ZonePtr).PeopelCountCalc_C = false;

                    helper_InternalThermalMassCalc = false;
                    helper_AirInfiltrationCalc = false;
                    helper_PeopleCountCalc = false;

                    helper_IncludeSystemSupply = false;

                    // Scenario 1: Only one unknown parameter to solve
                    // Scenario 1-1: To solve thermal mass
                    if (FlagHybridModel_TM && !FlagHybridModel_AI && !FlagHybridModel_PC) {
                        // Only measured temperature can be used to solve internal thermal mass
                        if (TemperatureSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = true;
                        }
                    }

                    // Scenario 1-2: To solve infiltration rate
                    if (!FlagHybridModel_TM && FlagHybridModel_AI && !FlagHybridModel_PC) {
                        // !!! Need to determine which variable to use if three variables are all available
                        // Temperature is measured
                        if (TemperatureSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_T = true;
                        }
                        // Humidity ratio is measured
                        if (HumidityRatioSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_H = true;
                        }
                        // CO2 concentration is measured
                        if (CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_C = true;
                        }
                    }

                    // Scenario 1-3: To solve people count
                    if (!FlagHybridModel_TM && !FlagHybridModel_AI && FlagHybridModel_PC) {
                        // !!! Need to determine which variable to use if three variables are all available
                        // Temperature is measured
                        if (TemperatureSchPtr > 0) {
                            HybridModelZone(ZonePtr).PeopelCountCalc_T = true;
                        }
                        // Humidity ratio is measured
                        if (HumidityRatioSchPtr > 0) {
                            HybridModelZone(ZonePtr).PeopelCountCalc_H = true;
                        }
                        // CO2 concentration is measured
                        if (CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).PeopelCountCalc_C = true;
                        }
                    }

                    // Scenario 2: Two unknown parameters to solve
                    // Scenario 2-1: To solve thermal mass and infiltration rate
                    if (FlagHybridModel_TM && FlagHybridModel_AI && !FlagHybridModel_PC) {
                        // Temperature and humidity ratio are measured
                        if (TemperatureSchPtr > 0 && HumidityRatioSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_H = true;
                        }
                        // Temperature and CO2 concentration are measured
                        if (TemperatureSchPtr > 0 && CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_C = true;
                        }
                    }

                    // Scenario 2-2: To solve thermal mass and people count
                    if (FlagHybridModel_TM && !FlagHybridModel_AI && FlagHybridModel_PC) {
                        // Temperature and humidity ratio are measured
                        if (TemperatureSchPtr > 0 && HumidityRatioSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_H = true;
                        }
                        // Temperature and CO2 concentration are measured
                        if (TemperatureSchPtr > 0 && CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_C = true;
                        }
                    }

                    // Scenario 2-3: To solve infiltration rate and people count
                    if (!FlagHybridModel_TM && FlagHybridModel_AI && FlagHybridModel_PC) {

                        // !!! Consider: which variable to use?

                        // Temperature and humidity ratio are measured
                        if (TemperatureSchPtr > 0 && HumidityRatioSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_T = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_T = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_H = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_H = true;
                        }
                        // Temperature and CO2 concentration are measured
                        if (TemperatureSchPtr > 0 && CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_T = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_T = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_C = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_C = true;
                        }
                        //  Humidity ratio and CO2 concentration are measured
                        if (HumidityRatioSchPtr > 0 && CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InfiltrationCalc_H = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_H = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_C = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_C = true;
                        }
                    }

                    // Scenario 3: Three unknown parameters to solve
                    if (FlagHybridModel_TM && FlagHybridModel_AI && FlagHybridModel_PC) {
                        if (TemperatureSchPtr > 0 && HumidityRatioSchPtr > 0 && CO2ConcentrationSchPtr > 0) {
                            HybridModelZone(ZonePtr).InternalThermalMassCalc_T = false;
                            HybridModelZone(ZonePtr).InfiltrationCalc_T = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_H = true;
                            HybridModelZone(ZonePtr).InfiltrationCalc_C = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_T = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_H = true;
                            HybridModelZone(ZonePtr).PeopelCountCalc_C = true;
                        }
                    }

                    // Summarise hybridmodel flags
                    if (HybridModelZone(ZonePtr).InternalThermalMassCalc_T) {
                        helper_InternalThermalMassCalc = true;
                    }
                    if (HybridModelZone(ZonePtr).InfiltrationCalc_T || HybridModelZone(ZonePtr).InfiltrationCalc_H ||
                        HybridModelZone(ZonePtr).InfiltrationCalc_C) {
                        helper_AirInfiltrationCalc = true;
                    }
                    if (HybridModelZone(ZonePtr).PeopelCountCalc_T || HybridModelZone(ZonePtr).PeopelCountCalc_H ||
                        HybridModelZone(ZonePtr).PeopelCountCalc_C) {
                        helper_PeopleCountCalc = true;
                    }

                    // Decide if system supply terms are sufficient to be included in the inverse solution
                    // For now, the inverse algorithms can only solve one unkown parameters at a time.
                    if ((SupplyAirTemperatureSchPtr > 0 && SupplyAirMassFlowRateSchPtr > 0) ||
                        (SupplyAirHumidityRatioSchPtr > 0 && SupplyAirMassFlowRateSchPtr > 0) ||
                        (SupplyAirCO2ConcentrationSchPtr > 0 && SupplyAirMassFlowRateSchPtr > 0)) {
                        HybridModelZone(ZonePtr).IncludeSystemSupplyParameters = true;
                        helper_IncludeSystemSupply = true;
                    }

                    // Need to update this section when writing the code to solve 2 or more variables at the same time
                    if (HybridModelZone(ZonePtr).InternalThermalMassCalc_T && HybridModelZone(ZonePtr).InfiltrationCalc_T) {
                        HybridModelZone(ZonePtr).InfiltrationCalc_T = false;
                        ShowWarningError(CurrentModuleObject + "=\"" + HybridModelZone(ZonePtr).Name + "\" invalid " + cAlphaFieldNames(3) + " and " +
                                         cAlphaFieldNames(4) + ".");
                        ShowContinueError("Field " + cAlphaFieldNames(3) + " and " + cAlphaFieldNames(4) + "\" cannot be both set to YES.");
                        ShowContinueError("Field " + cAlphaFieldNames(4) + "\" is changed to NO for the hybrid modeling simulations.");
                    }

                    // Flags showing Hybrid Modeling settings
                    // Need to update this if{}
                    FlagHybridModel = helper_InternalThermalMassCalc || helper_AirInfiltrationCalc || helper_PeopleCountCalc;

                    if (FlagHybridModel) {
                        // Here are the logic of getting inputs from idf.
                        // Need to update the section below
                        HybridModelZone(ZonePtr).ZoneMeasuredTemperatureSchedulePtr = GetScheduleIndex(cAlphaArgs(6));
                        HybridModelZone(ZonePtr).ZoneMeasuredHumidityRatioSchedulePtr = GetScheduleIndex(cAlphaArgs(7));
                        HybridModelZone(ZonePtr).ZoneMeasuredCO2ConcentrationSchedulePtr = GetScheduleIndex(cAlphaArgs(8));
                        HybridModelZone(ZonePtr).ZonePeopleActivityLevelSchedulePtr = GetScheduleIndex(cAlphaArgs(9));
                        HybridModelZone(ZonePtr).ZonePeopleSensibleFractionSchedulePtr = GetScheduleIndex(cAlphaArgs(10));
                        HybridModelZone(ZonePtr).ZonePeopleRadiationFractionSchedulePtr = GetScheduleIndex(cAlphaArgs(11));
                        HybridModelZone(ZonePtr).ZonePeopleCO2GenRateSchedulePtr = GetScheduleIndex(cAlphaArgs(12));

                        if (helper_IncludeSystemSupply) {
                            HybridModelZone(ZonePtr).ZoneSupplyAirTemperatureSchedulePtr = GetScheduleIndex(cAlphaArgs(13));
                            HybridModelZone(ZonePtr).ZoneSupplyAirMassFlowRateSchedulePtr = GetScheduleIndex(cAlphaArgs(14));
                            HybridModelZone(ZonePtr).ZoneSupplyAirHumidityRatioSchedulePtr = GetScheduleIndex(cAlphaArgs(15));
                            HybridModelZone(ZonePtr).ZoneSupplyAirCO2ConcentrationSchedulePtr = GetScheduleIndex(cAlphaArgs(16));
                        }

                        HybridModelZone(ZonePtr).ZoneMeasuredTemperatureStartMonth = rNumericArgs(1);
                        HybridModelZone(ZonePtr).ZoneMeasuredTemperatureStartDate = rNumericArgs(2);
                        HybridModelZone(ZonePtr).ZoneMeasuredTemperatureEndMonth = rNumericArgs(3);
                        HybridModelZone(ZonePtr).ZoneMeasuredTemperatureEndDate = rNumericArgs(4);

                        // prepare start and end date for Hybrid Modeling
                        {
                            int HMDayArr[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

                            HybridModelStartMonth = HybridModelZone(ZonePtr).ZoneMeasuredTemperatureStartMonth;
                            HybridModelStartDate = HybridModelZone(ZonePtr).ZoneMeasuredTemperatureStartDate;
                            HybridModelEndMonth = HybridModelZone(ZonePtr).ZoneMeasuredTemperatureEndMonth;
                            HybridModelEndDate = HybridModelZone(ZonePtr).ZoneMeasuredTemperatureEndDate;

                            if (HybridModelStartMonth >= 1 && HybridModelStartMonth <= 12) {
                                HMStartDay = HMDayArr[HybridModelStartMonth - 1];
                            } else {
                                HMStartDay = 0;
                            }

                            if (HybridModelEndMonth >= 1 && HybridModelEndMonth <= 12) {
                                HMEndDay = HMDayArr[HybridModelEndMonth - 1];
                            } else {
                                HMEndDay = 0;
                            }

                            HybridModelZone(ZonePtr).HybridStartDayOfYear = HMStartDay + HybridModelStartDate;
                            HybridModelZone(ZonePtr).HybridEndDayOfYear = HMEndDay + HybridModelEndDate;
                        }
                    }
                } else {
                    ShowSevereError(CurrentModuleObject + "=\"" + cAlphaArgs(1) + "\" invalid " + cAlphaFieldNames(2) + "=\"" + cAlphaArgs(2) +
                                    "\" not found.");
                    ErrorsFound = true;
                }
            }

            // ZoneAirMassFlowConservation should not be activated during the Hybrid Modeling infiltration calculations
            if (HybridModelZone(ZonePtr).InfiltrationCalc_T && ZoneAirMassFlow.EnforceZoneMassBalance) {
                ZoneAirMassFlow.EnforceZoneMassBalance = false;
                ShowWarningError("ZoneAirMassFlowConservation is deactivated when Hybrid Modeling is performed.");
            }

            // RoomAirModelType should be Mixing if Hybrid Modeling is performed for the zone
            if (FlagHybridModel) {
                for (ZonePtr = 1; ZonePtr <= NumOfZones; ZonePtr++) {
                    if ((HybridModelZone(ZonePtr).InternalThermalMassCalc_T || HybridModelZone(ZonePtr).InfiltrationCalc_T) &&
                        (AirModel(ZonePtr).AirModelType != RoomAirModel_Mixing)) {
                        AirModel(ZonePtr).AirModelType = RoomAirModel_Mixing;
                        ShowWarningError("Room Air Model Type should be Mixing if Hybrid Modeling is performed for the zone.");
                    }
                    if (helper_AirInfiltrationCalc) {
                        SetupOutputVariable("Zone Infiltration Hybrid Model Air Change Rate",
                                            OutputProcessor::Unit::ach,
                                            Zone(ZonePtr).InfilOAAirChangeRateHM,
                                            "Zone",
                                            "Average",
                                            Zone(ZonePtr).Name);
                        SetupOutputVariable("Zone Infiltration Hybrid Model Mass Flow Rate",
                                            OutputProcessor::Unit::kg_s,
                                            Zone(ZonePtr).MCPIHM,
                                            "Zone",
                                            "Average",
                                            Zone(ZonePtr).Name);
                    }
                    if (helper_PeopleCountCalc) {
                        SetupOutputVariable("Zone Hybrid Model People Count",
                                            OutputProcessor::Unit::None,
                                            Zone(ZonePtr).NumOccHM,
                                            "Zone",
                                            "Average",
                                            Zone(ZonePtr).Name);
                    }
                }
            }

            if (ErrorsFound) {
                ShowFatalError("Errors getting Hybrid Model input data. Preceding condition(s) cause termination.");
            }
        }
    }

    // Needed for unit tests, should not be normally called.
    void clear_state()
    {

        FlagHybridModel = false;
        NumOfHybridModelZones = 0;
    }

} // namespace HybridModel

} // namespace EnergyPlus
