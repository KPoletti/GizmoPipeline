import h5py
import numpy as np
import glob
from pytreegrav import ColumnDensity as ColumnDensity
import matplotlib.pyplot as plt


def getGasData(hdf5_file):
    """
    Reads in a hdf5 file and extracts the data into a dictionary.

    Args:
        hdf5_file (str): The path to the hdf5 file.
    output:
        data (dict): A dictionary containing the data from the hdf5 file.
    """
    file = h5py.File(hdf5_file, "r")
    data = {k: file["PartType0"][k][:] for k in file["PartType0"]}
    attr = {a: file["Header"].attrs[a] for a in file["Header"].attrs}
    return data, attr


def removeCellsOutsideSphere(data: dict, radius: float) -> np.ndarray:
    """
    Removes cells outside a sphere of a given radius.
    The sphere is centered at the origin of the coordinate system.

    Args:
        data (dict): data dictionary containing the coordinates of the cells
        radius (float): radius of the sphere to keep cells inside
    Returns:
        dict: data dictionary with cells outside the sphere removed
    """
    # center the origin
    cell_location = data["Coordinates"] - data["Coordinates"].max() / 2
    # calculate the distance from the center of the sphere
    distance = np.sqrt(np.sum(cell_location**2, axis=1))
    # keep only the cells inside the sphere
    index = np.where(distance <= radius)[0]
    # return the data dictionary with cells inside the sphere
    return {k: v[index] for k, v in data.items()}


def getTemperature(data: dict) -> np.ndarray:
    """
    Extracts the temperature from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The temperature array.
    """
    gamma = 1.200
    mean_molecular_weight = 1.22 * 1.67e-27  # kg
    k_boltzmann = 1.38064852e-23  # J/K
    energy = data["InternalEnergy"]
    return data["Temperature"]  # [code temperature] = [code energy]/[code mass]
    # return mean_molecular_weight * (gamma - 1) * energy / k_boltzmann


def getColumnDensity(data: dict, rays: int) -> np.ndarray:
    """
    Extracts the column density from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The column density array.
    """
    # Get the position and density arrays
    pos = data["Coordinates"]
    po_center = np.median(pos, axis=0)
    # Shift the position to the center of the box
    pos -= po_center
    mass = data["Masses"]
    density = data["Density"]
    volume = mass / density  # code units
    l_eff = (volume * (3 / 4 * np.pi)) ** (1 / 3)  # code units

    NH = ColumnDensity(
        pos, mass, l_eff, parallel=True, randomize_rays=False, rays=rays
    )  # [code mass]/[code length]^2
    return NH


def getRadiationDensity(data: dict, attr: dict) -> np.ndarray:
    """
    Extracts the radiation from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The UV radiation array.
    """
    # calculate the wavelength of the radiation field
    rad_min = attr["Radiation_RHD_Min_Bin_Freq_in_eV"]
    rad_max = attr["Radiation_RHD_Max_Bin_Freq_in_eV"]

    # select the photon energy as the radiation field
    rad = data["PhotonEnergy"]
    volume = np.zeros((rad.shape[0], 1))  # initialize the volume array
    # calculate the volume of the cell
    volume[:, 0] = data["Masses"] / data["Density"]  # code units
    # [code energy]/[code volume] = [code mass][code velocity]^2/[code length]^3
    return rad / volume


import h5py
import numpy as np
import glob
from pytreegrav import ColumnDensity as ColumnDensity
import matplotlib.pyplot as plt


def getGasData(hdf5_file):
    """
    Reads in a hdf5 file and extracts the data into a dictionary.

    Args:
        hdf5_file (str): The path to the hdf5 file.
    output:
        data (dict): A dictionary containing the data from the hdf5 file.
    """
    file = h5py.File(hdf5_file, "r")
    data = {k: file["PartType0"][k][:] for k in file["PartType0"]}
    attr = {a: file["Header"].attrs[a] for a in file["Header"].attrs}
    return data, attr


def removeCellsOutsideSphere(data: dict, radius: float) -> np.ndarray:
    """
    Removes cells outside a sphere of a given radius.
    The sphere is centered at the origin of the coordinate system.

    Args:
        data (dict): data dictionary containing the coordinates of the cells
        radius (float): radius of the sphere to keep cells inside
    Returns:
        dict: data dictionary with cells outside the sphere removed
    """
    # center the origin
    cell_location = data["Coordinates"] - data["Coordinates"].max() / 2
    # calculate the distance from the center of the sphere
    distance = np.sqrt(np.sum(cell_location**2, axis=1))
    # keep only the cells inside the sphere
    index = np.where(distance <= radius)[0]
    # return the data dictionary with cells inside the sphere
    return {k: v[index] for k, v in data.items()}


def getTemperature(data: dict) -> np.ndarray:
    """
    Extracts the temperature from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The temperature array.
    """
    gamma = 1.200
    mean_molecular_weight = 1.22 * 1.67e-27  # kg
    k_boltzmann = 1.38064852e-23  # J/K
    energy = data["InternalEnergy"]
    return data["Temperature"]  # [code temperature] = [code energy]/[code mass]
    # return mean_molecular_weight * (gamma - 1) * energy / k_boltzmann


def getColumnDensity(data: dict, rays: int) -> np.ndarray:
    """
    Extracts the column density from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The column density array.
    """
    # Get the position and density arrays
    pos = data["Coordinates"]
    po_center = np.median(pos, axis=0)
    # Shift the position to the center of the box
    pos -= po_center
    mass = data["Masses"]
    density = data["Density"]
    volume = mass / density  # code units
    l_eff = (volume * (3 / 4 * np.pi)) ** (1 / 3)  # code units

    NH = ColumnDensity(
        pos, mass, l_eff, parallel=True, randomize_rays=False, rays=rays
    )  # [code mass]/[code length]^2
    return NH


def getRadiationDensity(data: dict, attr: dict) -> np.ndarray:
    """
    Extracts the radiation from the data dictionary.

    Args:
        data (dict): The data dictionary.

    Returns:
        np.ndarray: The UV radiation array.
    """
    # calculate the wavelength of the radiation field
    rad_min = attr["Radiation_RHD_Min_Bin_Freq_in_eV"]
    rad_max = attr["Radiation_RHD_Max_Bin_Freq_in_eV"]

    # select the photon energy as the radiation field
    rad = data["PhotonEnergy"]
    volume = np.zeros((rad.shape[0], 1))  # initialize the volume array
    # calculate the volume of the cell
    volume[:, 0] = data["Masses"] / data["Density"]  # code units
    # [code energy]/[code volume] = [code mass][code velocity]^2/[code length]^3
    return rad / volume
