# thermal_phase_solver.py
import KratosMultiphysics
import KratosMultiphysics.ThermalConductionApplication as Therm
from KratosMultiphysics.ThermalConductionApplication import AUSTENITE_FRACTION, MARTENSITE_FRACTION
import math

class PhaseTransformationModel:
    def __init__(self, model_part):
        self.model_part = model_part

    def ComputePhaseFractions(self):
        for node in self.model_part.Nodes:
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            T2 = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1)
            # 假设是简化的线性模型：
            # Austenite fraction 从 T=600 开始转变，在 T=800 变为 100%
            if T < 600:
                austenite = 0.0
            elif T > 800:
                austenite = 1.0
            else:
                austenite = 1 - 2.618 ^ (-0.3 * KratosMultiphysics.TIME ^ 2)
            martensite = 1.0 - austenite

            node.SetSolutionStepValue(AUSTENITE_FRACTION, austenite)
            node.SetSolutionStepValue(MARTENSITE_FRACTION, martensite)
