import streamlit as st
from functions import compute_Z_Hall_Yarborough, compute_Z_with_DAK
from DPR import compute_Z_DPR

def main():
    st.title("Z-Factor Calculator")

    fluid_type = st.selectbox("Fluid Type", ["NATURAL GAS", "CONDENSATE"])
    sp_gravity = st.number_input("Specific Gravity", min_value=0.0, format="%.4f")
    temperature = st.number_input("Temperature (Rankine)", min_value=0.0, format="%.2f")
    pressure = st.number_input("Pressure (psia)", min_value=0.0, format="%.2f")
    method = st.selectbox("Calculation Method", ["Hall-Yarborough", "Dranchuk-Abou-Kaseem", "Dranchuk-Purvis-Robinson"])

    if st.button("Calculate"):
        try:
            if method == "Hall-Yarborough":
                z_factor, Pc, Tc = compute_Z_Hall_Yarborough(sp_gravity, temperature, pressure, fluid_type)
            elif method == "Dranchuk-Abou-Kaseem":
                z_factor, Pc, Tc = compute_Z_with_DAK(sp_gravity, temperature, pressure, fluid_type)
            else:
                z_factor, Pc, Tc = compute_Z_DPR(sp_gravity, temperature, pressure, fluid_type)

            st.success(f"The pseudo-critical temperature of the fluid is: {Tc}")
            st.success(f"The pseudo-critical pressure of the fluid is: {Pc}")
            st.success(f"The Z-Factor is: {z_factor}")
        except ValueError as e:
            st.error(f"Error: {e}")
        except Exception as e:
            st.error(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
