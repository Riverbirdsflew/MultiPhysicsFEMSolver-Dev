import KratosMultiphysics as KM
import KratosMultiphysics.StructuralAnalysisApplication as KSA


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    model_part = Model.GetModelPart(settings["Parameters"]["model_part_name"].GetString())
    return KSA.DistributeLoadOnSurfaceProcess(model_part, settings["Parameters"])
