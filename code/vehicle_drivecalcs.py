import math

MASS_OF_VEHICLE_KEY = "mass"
DRAG_COEFFICIENT_KEY = "drag_coefficient"
FRONTAL_AREA_OF_VEHICLE_KEY = "area_vehicle"
ROLLING_RESISTANCE_KEY = "rolling_resistance"
ACCESSORY_POWER_DEMAND_KEY = "accessory_power_demand"
ENGINE_EFFICIENCY_ALPHA_CONSTANT_KEY = "eta_alpha_constant"
ENGINE_EFFICIENCY_WORK_CONSTANT_1_KEY = "eta_w1_constant"
ENGINE_EFFICIENCY_WORK_CONSTANT_2_KEY = "eta_w2_constant"

GASOLINE_ENERGY_DENSITY_KEY = "gasoline_energy_density"
AIR_VELOCITY_KEY = "air_velocity"
AIR_DENSITY_KEY = "air_density"
ROAD_ANGLE_KEY = "road_angle"
GRAVITATIONAL_CONSTANT_KEY = "gravitational_constant"

TIME_COLNAME = "time_s"
VELOCITY_VEHICLE_MPH_COLNAME = "v_mph"
VELOCITY_VEHICLE_M_PER_S_COLNAME = "v_veh_m_per_s"
ACCELERATION_M_PER_S_COLNAME = "accel_m_per_s"
DELTA_DISTANCE_KM_COLNAME = "delta_distance_km"
TOTAL_DISTANCE_KM_COLNAME = "total_distance_km"
POWER_DEMAND_BRAKING_COLNAME = "brake_force_n"
POWER_DEMAND_DRAG_COLNAME = "power_demand_drag_kw"
POWER_DEMAND_ROLL_COLNAME = "power_demand_rolling_resistance_kw"
POWER_DEMAND_ACCEL_COLNAME = "power_demand_acceleration_kw"
POWER_DEMAND_TOTAL_COLNAME = "power_demand_vehicle_kw"
POWER_DEMAND_ACCESSORIES_COLNAME = "power_demand_accessories_kw"
TOTAL_VEHICLE_ENERGY_COLNAME = "total_vehicle_energy_mj"
AVG_ENERGY_VEHICLE_PER_VKT_COLNAME = "average_e_veh_over_vkt_mj_per_km"
DRIVETRAIN_EFFICIENCY_COLNAME = "eta"
INSTANTANEOUS_FUEL_REQT_COLNAME = "instant_fuel_req_kw"
TOTAL_FUEL_ENERGY_COLNAME = "total_fuel_energy_mj"
FUEL_CONSUMPTION_LITERS_COLNAME = "fuel_consumption_liters"
FUEL_EFFICIENCY_KM_PER_L_COLNAME = "fuel_efficiency_km_per_l"
AVG_FUEL_MJ_PER_KM_COLNAME = "avg_fuel_mj_per_km"

TOTAL_DISTANCE_TRAVELED_KEY = "total_distance_traveled"
MAX_VEHICLE_POWER_DEMAND_KEY = "max_vehicle_power_demand"
AVG_VEHICLE_POWER_DEMAND_KEY = "avg_vehicle_power_demand"
TOTAL_VEHICLE_ENERGY_DEMAND_KEY = "total_vehicle_energy_demand"
TOTAL_FUEL_ENERGY_DEMAND_KEY = "total_fuel_energy_demand"
TOTAL_FUEL_CONSUMED_KEY = "total_fuel_consumed"
AVG_FUEL_ECONOMY_LM_PER_L_KEY = "avg_fuel_economy_km_per_l"
AVG_FUEL_ECONOMY_MPG_KEY = "avg_fuel_economy_mpg"

def get_instantaneous_drivetrain_efficiency(alpha, work_vehicle, w_1, w_2):
    # what does w_1 and w_2 represent?
    
    inner_term_one = math.exp( - (work_vehicle / w_1))
         
    inner_term_two = math.exp( - (work_vehicle / w_2))

    inner_term = inner_term_one - inner_term_two
    
    return alpha * inner_term

def get_m_per_s_from_mph(mph):
    
    return mph * 1609.344 / 3600

def get_denominator(unit):
    
    valid_units = ["kW", "W"]
    if unit not in valid_units:
        raise ValueError("Unit to return ({}) not one of valid options: {}".format(unit, valid_units))
    
    denominator = 1. if unit == "W" else 1000.
    
    return denominator

def get_power_demand_drag(drag_coefficient, frontal_area_of_vehicle, air_density, vehicle_velocity_m_per_s, air_velocity, unit = "W"):
    """
    Gets power demand of vehicle due to vehicle drag
    
    get_power_demand_drag(0.33, 2.2, 1.149, get_m_per_s_from_mph(3.), 0, "kW") = 0.001006
    """
    
    valid_units = ["kW", "W"]
    if unit not in valid_units:
        raise ValueError("Unit to return ({}) not one of valid options: {}".format(unit, valid_units))
    
    denominator = 1. if unit == "W" else 1000.
    return (
        drag_coefficient \
        * frontal_area_of_vehicle \
        * (air_density / 2) \
        * (vehicle_velocity_m_per_s - air_velocity)**2
        ) \
        * vehicle_velocity_m_per_s  \
    / denominator

def get_power_demand_roll(rolling_resistance_coefficient, vehicle_mass_kg, gravitational_constant, vehicle_velocity_m_per_s, unit = "W"):
    """
    Gets power demand of vehicle due to rolling resistance 
    
    get_power_demand_roll(0.009, 3000., GRAVITATIONAL_CONSTANT, get_m_per_s_from_mph(8.6), "kW") = 1.0183
    """
    valid_units = ["kW", "W"]
    if unit not in valid_units:
        raise ValueError("Unit to return ({}) not one of valid options: {}".format(unit, valid_units))
    
    denominator = 1. if unit == "W" else 1000.
    
    return (
        vehicle_mass_kg \
        * rolling_resistance_coefficient \
        * gravitational_constant \
        * vehicle_velocity_m_per_s
    ) / denominator

def get_power_demand_acceleration(vehicle_mass_kg, vehicle_velocity_m_per_s, vehicle_acceleration_m_per_s, unit = "W"):
    
    denominator = get_denominator(unit)
    return (
        vehicle_mass_kg \
        * vehicle_acceleration_m_per_s \
        * vehicle_velocity_m_per_s
    ) / denominator

def get_braking_force_n(
    drag_coefficient, 
    rolling_resistance_coefficient,
    vehicle_mass_kg, 
    frontal_area_of_vehicle, 
    gravitational_constant,
    air_density, 
    air_velocity, 
    velocity_vehicle_m_per_s, 
    acceleration_vehicle_m_per_s,
):
    """
    get_braking_force_n(
        drag_coefficient=0.33,
        rolling_resistance_coefficient=0.009,
        vehicle_mass_kg=3000,
        frontal_area_of_vehicle=2.2,
        gravitational_constant=GRAVITATIONAL_CONSTANT,
        air_density=1.149,
        air_velocity=0,
        velocity_vehicle_m_per_s=9.61,
        acceleration_vehicle_m_per_s=-0.268224
    ) = 578.3208
    """
    braking_term_one = (
        drag_coefficient \
        * frontal_area_of_vehicle \
        / 2.
    )
    
    braking_term_two = (
        air_density \
        * (velocity_vehicle_m_per_s - air_velocity) ** 2
    )
    
    braking_term_three = (
        vehicle_mass_kg \
        * (rolling_resistance_coefficient * gravitational_constant + acceleration_vehicle_m_per_s)
    )
    
    return max(0, (braking_term_one * braking_term_two) - braking_term_three)

def create_cycle_df(df, vehicle_config, constants_config):
    """
    Controller function to create all columns needed for a cycle given time and velocity columns
    """
    cycle_df = df.copy()
    
    cycle_df[VELOCITY_VEHICLE_M_PER_S_COLNAME] = cycle_df.apply(
        lambda x: get_m_per_s_from_mph(x[VELOCITY_VEHICLE_MPH_COLNAME]),
        axis = 1
    )
    
    acceleration_array = []
    delta_distance_km_array = []
    total_distance_km_array = []
    
    for idx, row in cycle_df.iterrows():

        if idx == 0:
            acceleration = 0
            delta_distance_km = 0
            total_distance_km = 0

        else:
            acceleration = cycle_df.iloc[idx][VELOCITY_VEHICLE_M_PER_S_COLNAME] - cycle_df.iloc[idx - 1][VELOCITY_VEHICLE_M_PER_S_COLNAME] 
            delta_distance_km = ((
                cycle_df.iloc[idx][VELOCITY_VEHICLE_M_PER_S_COLNAME] \
                + cycle_df.iloc[idx - 1][VELOCITY_VEHICLE_M_PER_S_COLNAME] 
                ) / 2 \
                * (cycle_df.iloc[idx][TIME_COLNAME] - cycle_df.iloc[idx - 1][TIME_COLNAME])
                ) / 1000

        acceleration_array.append(acceleration)
        delta_distance_km_array.append(delta_distance_km)

    cycle_df[ACCELERATION_M_PER_S_COLNAME] = acceleration_array
    cycle_df[DELTA_DISTANCE_KM_COLNAME] = delta_distance_km_array
    
    for idx, row in cycle_df.iterrows():
        
        total_distance_km_array.append(sum(cycle_df[0:idx][DELTA_DISTANCE_KM_COLNAME]) + row[DELTA_DISTANCE_KM_COLNAME])
        
    cycle_df[TOTAL_DISTANCE_KM_COLNAME] = total_distance_km_array
    
    cycle_df[POWER_DEMAND_BRAKING_COLNAME] = cycle_df.apply(
        lambda x: get_braking_force_n(
            drag_coefficient= vehicle_config[DRAG_COEFFICIENT_KEY],
            rolling_resistance_coefficient= vehicle_config[ROLLING_RESISTANCE_KEY],
            vehicle_mass_kg= vehicle_config[MASS_OF_VEHICLE_KEY],
            frontal_area_of_vehicle= vehicle_config[FRONTAL_AREA_OF_VEHICLE_KEY],
            gravitational_constant= constants_config[GRAVITATIONAL_CONSTANT_KEY],
            air_density=constants_config[AIR_DENSITY_KEY],
            air_velocity=constants_config[AIR_VELOCITY_KEY],
            velocity_vehicle_m_per_s=x[VELOCITY_VEHICLE_M_PER_S_COLNAME],
            acceleration_vehicle_m_per_s=x[ACCELERATION_M_PER_S_COLNAME]
        ),
        axis = 1
    )

    cycle_df[POWER_DEMAND_DRAG_COLNAME] = cycle_df.apply(
        lambda x: 0 if x[POWER_DEMAND_BRAKING_COLNAME] > 0 else get_power_demand_drag(
            drag_coefficient = vehicle_config[DRAG_COEFFICIENT_KEY],
            frontal_area_of_vehicle = vehicle_config[FRONTAL_AREA_OF_VEHICLE_KEY],
            air_density=constants_config[AIR_DENSITY_KEY],
            vehicle_velocity_m_per_s= x[VELOCITY_VEHICLE_M_PER_S_COLNAME],
            air_velocity= constants_config[AIR_VELOCITY_KEY],
            unit = "kW"
        ),
        axis=1
    )

    cycle_df[POWER_DEMAND_ROLL_COLNAME] = cycle_df.apply(
        lambda x: 0 if x[POWER_DEMAND_BRAKING_COLNAME] > 0 else get_power_demand_roll(
            rolling_resistance_coefficient= vehicle_config[ROLLING_RESISTANCE_KEY],
            vehicle_mass_kg= vehicle_config[MASS_OF_VEHICLE_KEY],
            gravitational_constant= constants_config[GRAVITATIONAL_CONSTANT_KEY],
            vehicle_velocity_m_per_s= x[VELOCITY_VEHICLE_M_PER_S_COLNAME],
            unit = "kW"
        ),
        axis = 1
    )

    cycle_df[POWER_DEMAND_ACCEL_COLNAME] = cycle_df.apply(
        lambda x: 0 if x[POWER_DEMAND_BRAKING_COLNAME] > 0 else get_power_demand_acceleration(
            vehicle_mass_kg= vehicle_config[MASS_OF_VEHICLE_KEY],
            vehicle_velocity_m_per_s= x[VELOCITY_VEHICLE_M_PER_S_COLNAME],
            vehicle_acceleration_m_per_s= x[ACCELERATION_M_PER_S_COLNAME],
            unit = "kW"
        ),
        axis = 1
    )

    cycle_df[POWER_DEMAND_ACCESSORIES_COLNAME] = vehicle_config[ACCESSORY_POWER_DEMAND_KEY]

    cycle_df[POWER_DEMAND_TOTAL_COLNAME] = cycle_df[POWER_DEMAND_ACCEL_COLNAME] \
        + cycle_df[POWER_DEMAND_ACCESSORIES_COLNAME] \
        + cycle_df[POWER_DEMAND_DRAG_COLNAME] \
        + cycle_df[POWER_DEMAND_ROLL_COLNAME]

    cumulative_energy_mj = 0
    cumulative_energy_array = []

    for idx, row in cycle_df.iterrows():
        if idx == 0:
            cumulative_energy_array.append(0)
            continue

        time_diff = cycle_df.iloc[idx][TIME_COLNAME] - cycle_df.iloc[idx - 1][TIME_COLNAME] 
        power_demand_sum = cycle_df.iloc[idx][POWER_DEMAND_TOTAL_COLNAME] + cycle_df.iloc[idx - 1][POWER_DEMAND_TOTAL_COLNAME]
        addtl_energy = time_diff * power_demand_sum / 2000

        cumulative_energy_mj += addtl_energy
        cumulative_energy_array.append(cumulative_energy_mj)

    cycle_df[TOTAL_VEHICLE_ENERGY_COLNAME] = cumulative_energy_array

    cycle_df[AVG_ENERGY_VEHICLE_PER_VKT_COLNAME] = cycle_df.apply(
        lambda x: x[TOTAL_VEHICLE_ENERGY_COLNAME] / x[TOTAL_DISTANCE_KM_COLNAME] if x[TOTAL_DISTANCE_KM_COLNAME] > 0 else 0,
        axis = 1
    )

    cycle_df[DRIVETRAIN_EFFICIENCY_COLNAME] = cycle_df.apply(
        lambda x: get_instantaneous_drivetrain_efficiency(
            alpha=vehicle_config[ENGINE_EFFICIENCY_ALPHA_CONSTANT_KEY],
            work_vehicle= x[POWER_DEMAND_TOTAL_COLNAME],
            w_1=vehicle_config[ENGINE_EFFICIENCY_WORK_CONSTANT_1_KEY],
            w_2=vehicle_config[ENGINE_EFFICIENCY_WORK_CONSTANT_2_KEY]
        ),
        axis = 1
    )

    cycle_df[INSTANTANEOUS_FUEL_REQT_COLNAME] = cycle_df.apply(
        lambda x: x[POWER_DEMAND_TOTAL_COLNAME] / x[DRIVETRAIN_EFFICIENCY_COLNAME],
        axis = 1
    )
    
    cumulative_fuel_energy_array = []
    cumulative_fuel_energy_mj = 0

    for idx, row in cycle_df.iterrows():
        if idx == 0:
            cumulative_fuel_energy_array.append(0)
            continue

        time_diff = cycle_df.iloc[idx][TIME_COLNAME] - cycle_df.iloc[idx - 1][TIME_COLNAME] 
        fuel_req = (cycle_df.iloc[idx][INSTANTANEOUS_FUEL_REQT_COLNAME] + cycle_df.iloc[idx - 1][INSTANTANEOUS_FUEL_REQT_COLNAME]) / 2.

        cumulative_fuel_energy_mj += (time_diff * fuel_req) / 1000

        cumulative_fuel_energy_array.append(cumulative_fuel_energy_mj)

    cycle_df[TOTAL_FUEL_ENERGY_COLNAME] = cumulative_fuel_energy_array
    
    cycle_df[FUEL_CONSUMPTION_LITERS_COLNAME] = cycle_df.apply(
        lambda x: x[TOTAL_FUEL_ENERGY_COLNAME] / constants_config[GASOLINE_ENERGY_DENSITY_KEY],
        axis = 1
    )

    cycle_df[FUEL_EFFICIENCY_KM_PER_L_COLNAME] = cycle_df.apply(
        lambda x: None if x[FUEL_CONSUMPTION_LITERS_COLNAME] == 0 \
            else x[TOTAL_DISTANCE_KM_COLNAME] / x[FUEL_CONSUMPTION_LITERS_COLNAME],
        axis = 1
    )

    cycle_df[AVG_FUEL_MJ_PER_KM_COLNAME] = cycle_df.apply(
        lambda x: None if x[FUEL_EFFICIENCY_KM_PER_L_COLNAME] in [None, 0] \
            else x[TOTAL_FUEL_ENERGY_COLNAME] / x[TOTAL_DISTANCE_KM_COLNAME],
        axis = 1
    )
    
    return cycle_df

def find_aggregate_stats(cycle_df):
    
    total_distance_traveled = max(cycle_df[TOTAL_DISTANCE_KM_COLNAME])
    max_vehicle_power_demand = max(cycle_df[POWER_DEMAND_TOTAL_COLNAME])
    avg_vehicle_power_demand = max(cycle_df[TOTAL_VEHICLE_ENERGY_COLNAME]) / max(cycle_df[TIME_COLNAME]) * 1000
    total_vehicle_energy_demand = max(cycle_df[TOTAL_VEHICLE_ENERGY_COLNAME])
    total_fuel_energy_demand = max(cycle_df[TOTAL_FUEL_ENERGY_COLNAME])
    total_fuel_consumed = max(cycle_df[FUEL_CONSUMPTION_LITERS_COLNAME])
    avg_fuel_economy_km_per_l = total_distance_traveled / max(cycle_df[FUEL_CONSUMPTION_LITERS_COLNAME])
    avg_fuel_economy_mpg = avg_fuel_economy_km_per_l * 2.35215
    
    return {
        TOTAL_DISTANCE_TRAVELED_KEY : total_distance_traveled,
        MAX_VEHICLE_POWER_DEMAND_KEY : max_vehicle_power_demand,
        AVG_VEHICLE_POWER_DEMAND_KEY : avg_vehicle_power_demand,
        TOTAL_VEHICLE_ENERGY_DEMAND_KEY : total_vehicle_energy_demand,
        TOTAL_FUEL_ENERGY_DEMAND_KEY : total_fuel_energy_demand,
        TOTAL_FUEL_CONSUMED_KEY : total_fuel_consumed,
        AVG_FUEL_ECONOMY_LM_PER_L_KEY : avg_fuel_economy_km_per_l,
        AVG_FUEL_ECONOMY_MPG_KEY : avg_fuel_economy_mpg
    }