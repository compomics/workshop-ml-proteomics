import numpy as np
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
from pyteomics.usi import USI, proxi
from pyteomics.mass.unimod import Unimod
from ms2pip.single_prediction import SinglePrediction

def reshape_attributes(attributes):
    return {item["name"]: item["value"] for item in attributes}


def get_usi_spectrum(usi, **spectrum_utils_kwargs):
    proxi_response = proxi(usi)
    proxi_attributes = reshape_attributes(proxi_response["attributes"])

    # Process spectrum and plot
    spectrum = (
        sus.MsmsSpectrum(
            identifier=usi,
            **spectrum_utils_kwargs,
            precursor_mz=proxi_attributes["isolation window target m/z"],
            precursor_charge=int(proxi_attributes["charge state"]),
            mz=proxi_response["m/z array"],
            intensity=proxi_response["intensity array"],
            peptide=proxi_attributes["unmodified peptide sequence"],
        )
        .filter_intensity(0.01, 50)
    )
    return spectrum


def get_theoretical_spectrum(peptide, modifications, b_ion_weight=1.0, ms2pip_instance=None):
    if not ms2pip_instance:
        ms2pip = SinglePrediction()
    else:
        ms2pip = ms2pip_instance

    charge = 2  # Does not really matter here

    mz, _, annotation = ms2pip.predict(peptide, modifications, charge, model="HCD2019")
    intensity = [b_ion_weight if ann[0] == "b" else 1.0 for ann in annotation]

    identifier = f"{peptide}/{charge}/{modifications}"
    precursor_mz = ms2pip.mod_info.calc_precursor_mz(peptide, modifications, charge)
    mod_dict = ms2pip._modifications_to_dict(modifications)
    sus_annotation = ms2pip._get_sus_annotation(mz, annotation)

    spectrum = sus.MsmsSpectrum(
        identifier,
        precursor_mz,
        charge,
        mz,
        intensity,
        annotation=sus_annotation,
        peptide=peptide,
        modifications=mod_dict,
    )
    return spectrum


def get_predicted_spectrum(peptide, modifications, charge, model="HCD2019", ms2pip_instance=None):
    if not ms2pip_instance:
        ms2pip = SinglePrediction()
    else:
        ms2pip = ms2pip_instance

    mz, intensity, annotation = ms2pip.predict(peptide, modifications, charge, model=model)

    identifier = f"{peptide}/{charge}/{modifications}"
    precursor_mz = ms2pip.mod_info.calc_precursor_mz(peptide, modifications, charge)
    mod_dict = ms2pip._modifications_to_dict(modifications)
    sus_annotation = ms2pip._get_sus_annotation(mz, annotation)

    spectrum = sus.MsmsSpectrum(
        identifier,
        precursor_mz,
        charge,
        mz,
        intensity,
        annotation=sus_annotation,
        peptide=peptide,
        modifications=mod_dict,
    )
    return spectrum


def _sus_annotation_to_index(ann, num_ions):
    """Convert spectrum_utils annotation to correct array index."""
    col_indices = {("b", 1): 0, ("y", 1): 1}
    if ann.ion_type == "b":
        return col_indices[ann.ion_type, ann.charge], int(ann.ion_index) - 1
    elif ann.ion_type == "y":
        return col_indices[ann.ion_type, ann.charge], num_ions - int(ann.ion_index)


def get_intensity_array(spectrum):
    n_ions = len(spectrum.peptide) - 1
    intensity = np.zeros((2, n_ions))
    x, y = np.vectorize(_sus_annotation_to_index)(
        spectrum.annotation[spectrum.annotation != None],
        n_ions
    )
    intensity_flat = np.vectorize(lambda x: x.calc_mz)(spectrum.annotation[spectrum.annotation != None])
    intensity[x, y] = intensity_flat
    return intensity