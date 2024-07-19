import streamlit as st
from functions import compute_Z_Hall_Yarborough, compute_Z_with_DAK
from DPR import compute_Z_DPR

def convert_to_rankine(value, unit):
    if unit == "Celsius":
        return (value + 273.15) * 1.8
    elif unit == "Fahrenheit":
        return value + 460
    elif unit == "Kelvin":
        return value * 1.8
    elif unit == "Rankine":
        return value
    else:
        raise ValueError("Unsupported temperature unit")


def convert_to_psia(value, unit):
    if unit == "psia":
        return value
    elif unit == "psig":
        return value + 14.7
    elif unit == "bar":
        return value * 14.5038
    elif unit == "kPa":
        return value * 0.145038
    else:
        raise ValueError("Unsupported pressure unit")


def main():
    st.title("Z-Factor Calculator")

    fluid_type = st.selectbox("Fluid Type", ["NATURAL GAS", "CONDENSATE"])
    sp_gravity = st.number_input("Specific Gravity", min_value=0.0, format="%.4f")

    temp_unit = st.selectbox("Temperature Unit", ["Celsius", "Fahrenheit", "Kelvin", "Rankine"])
    temperature = st.number_input(f"Temperature ({temp_unit})", min_value=0.0, format="%.2f")

    pressure_unit = st.selectbox("Pressure Unit", ["psia", "psig", "bar", "kPa"])
    pressure = st.number_input(f"Pressure ({pressure_unit})", min_value=0.0, format="%.2f")

    method = st.selectbox("Calculation Method", ["Hall-Yarborough", "Dranchuk-Abou-Kaseem", "Dranchuk-Purvis-Robinson"])

    if st.button("Calculate"):
        try:
            temp_in_rankine = convert_to_rankine(temperature, temp_unit)
            pressure_in_psia = convert_to_psia(pressure, pressure_unit)

            if method == "Hall-Yarborough":
                z_factor, Pc, Tc = compute_Z_Hall_Yarborough(sp_gravity, temp_in_rankine, pressure_in_psia, fluid_type)
            elif method == "Dranchuk-Abou-Kaseem":
                z_factor, Pc, Tc = compute_Z_with_DAK(sp_gravity, temp_in_rankine, pressure_in_psia, fluid_type)
            else:
                z_factor, Pc, Tc = compute_Z_DPR(sp_gravity, temp_in_rankine, pressure_in_psia, fluid_type)

            st.success(f"The pseudo-critical temperature of the fluid is: {Tc}")
            st.success(f"The pseudo-critical pressure of the fluid is: {Pc}")
            st.success(f"The Z-Factor is: {z_factor}")
        except ValueError as e:
            st.error(f"Error: {e}")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    main()
