import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
from pyteomics.usi import USI, proxi
from pyteomics.mass.unimod import Unimod
from ms2pip.single_prediction import SinglePrediction

def reshape_attributes(attributes):
    return {item["name"]: item["value"] for item in attributes}


def get_usi_spectrum(usi):
    proxi_response = proxi(usi)
    proxi_attributes = reshape_attributes(proxi_response["attributes"])

    # Process spectrum and plot
    spectrum = (
        sus.MsmsSpectrum(
            identifier=usi,
            precursor_mz=proxi_attributes["isolation window target m/z"],
            precursor_charge=int(proxi_attributes["charge state"]),
            mz=proxi_response["m/z array"],
            intensity=proxi_response["intensity array"],
            peptide=proxi_attributes["unmodified peptide sequence"],
            modifications=None
        )
        .filter_intensity(0.01, 50)
    )
    return spectrum


def get_theoretical_spectrum(peptide, modifications, charge, b_ion_weight=1.0, ms2pip_instance=None):
    if not ms2pip_instance:
        ms2pip = SinglePrediction()
    else:
        ms2pip = ms2pip_instance

    mz, _, annotation = ms2pip.predict(peptide, modifications, charge)
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


def get_predicted_spectrum(peptide, modifications, charge, ms2pip_instance=None):
    if not ms2pip_instance:
        ms2pip = SinglePrediction()
    else:
        ms2pip = ms2pip_instance
        
    mz, intensity, annotation = ms2pip.predict(peptide, modifications, charge)
    
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

