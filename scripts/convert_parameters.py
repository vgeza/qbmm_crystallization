#!/usr/bin/env python
# -*- coding: utf-8 -*-

# python -m pip install argparse pint PyYAML

"""
Script for converting PBM coefficients from literature to OpenFOAM format.
SI units are used for the converted set of parameters.
"""

import argparse
import pint
import yaml

VERSION = "1.1"
UREG = pint.UnitRegistry()


def default_value(name: str):
    """Get default value of a parameter, based on its unit. Returns zero."""
    unit = get_SI_unit(name)

    if unit is None:
        return 0

    return UREG.Quantity(f"0 {unit}")


def get_SI_unit(name: str) -> str:
    """Get SI unit of a parameter, based on its name."""
    if name == "R":
        return "J mol^-1 K^-1"

    if name.startswith("E"):
        return "J mol^-1"

    if name.startswith("alpha"):
        return "m^-1"

    if name.startswith("rho"):
        return "kg m^-3"

    return None


def convert_paramter(name, value):
    """Convert dimensionless numbers and simple dimensional quantities."""
    unit = get_SI_unit(name)

    if unit is None:
        return value

    return UREG.Quantity(value).to(unit)


def convert_paramter_list(names: list[str], values: list) -> dict:
    """Convert multiple parameters. For missing parameters, a default value of
    0 is assumed"""
    result = {}

    for p in names:
        value = values.get(p, default_value(p))
        result[p] = convert_paramter(p, value)

    return result


def convert_general_parameters(parameters_source: dict) -> dict:
    """Convert parameters in the 'general' section."""
    par = parameters_source.get("general", {})
    names = ["R", "rhoCrystal", "kv"]
    result = convert_paramter_list(names, par)

    # add comment as the first entry
    if "comment" in par:
        result = {"//": par["comment"]} | result

    return result


def convert_growth_parameters(parameters_source: dict) -> dict:
    """Convert parameters in the 'growth' section."""
    par = parameters_source.get("growth", {})
    names = [
        "g",
        "Eg",
        "alphaG",
        "Lg",
        "d",
        "Ed",
        "alphaD",
        "Ld",
        "relativeSupersaturation",
    ]
    result = convert_paramter_list(names, par)

    # add comment as the first entry
    if "comment" in par:
        result = {"//": par["comment"]} | result

    # TODO: convert kG and kD

    return result


def convert_nucleation_parameters(parameters_source: dict) -> dict:
    """Convert parameters in the 'nucleation' section."""
    par = parameters_source.get("nucleation", {})
    names = ["s", "Ms", "Es", "relativeSupersaturation"]
    result = convert_paramter_list(names, par)

    # add comment as the first entry
    if "comment" in par:
        result = {"//": par["comment"]} | result

    # TODO: convert kS

    return result


def write_OF_dictionaries(parameters: dict) -> None:
    """Write parameters as OpenFOAM dictionaries (generalDict, growthDict and
    nucleationDict). Header and '{}' are omitted since they are supposed to be
    #include'd in crystallizationProperties and populationBalanceProperties."""
    indent = " " * 8
    for section, par in parameters.items():
        filename = f"{section}Dict"
        print(f"Writing {filename}")
        with open(filename, "w") as f:
            for name, value in par.items():
                unit = get_SI_unit(name)
                if unit is None:
                    value_OF = value
                else:
                    value_OF = f"[{unit}] {value.magnitude}"
                f.write(f"{indent}{name} {value_OF};\n")


def main() -> None:
    """Read the input data, convert parameters and write results."""
    parser = argparse.ArgumentParser(
        prog="convert_parameters",
        description=(
            "Convert growth/nucleation YAML parameters into "
            "a validated OpenFOAM dictionary."
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {VERSION}"
    )
    parser.add_argument(
        "-s", "--source", default="parameters.yaml", help="Input YAML file"
    )
    parser.add_argument(
        "-d",
        "--destination",
        default="dict",
        help="Output OpenFOAM dictionary",
    )
    args = parser.parse_args()

    with open(args.source) as f:
        print(f"Reading {args.source}")
        parameters = yaml.safe_load(f)
    if not parameters:
        parameters = {}
    print(f"\n{parameters=}")

    parameters_converted = {}
    parameters_converted["general"] = convert_general_parameters(parameters)
    parameters_converted["growth"] = convert_growth_parameters(parameters)
    parameters_converted["nucleation"] = convert_nucleation_parameters(
        parameters
    )
    print(f"\n{parameters_converted=}\n")

    write_OF_dictionaries(parameters_converted)
    print("Done!")


if __name__ == "__main__":
    main()
