import numpy as np
import matplotlib.pyplot as plt
import math as math
import random as rnd
import itertools
from ortools.sat.python import cp_model
from datetime import datetime, timedelta
from math import log, exp
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import copy
import matplotlib.patches as patches




def time_formatting(time):
  #"""Takes the time information in form of tuple (year, month, day, hour, minute) and converts it to a datetime object"""
  delivery_year, delivery_month, delivery_day, delivery_hour, delivery_minute = time
  return(datetime(delivery_year, delivery_month, delivery_day, delivery_hour, delivery_minute, 0))


class Person:
  def __init__(self, id, time_constraints, mass, type):
    self.id = id
    self.time_constraints = time_constraints
    self.mass = mass
    self.accumulation_length = 150
    self.PET_length = 30
    self.treatment = "FDG"

#keeps track of 1 delivery (delivery_time, activity)
class Delivery:
  def __init__(self, time_delivery, delivery_activity):
    self.starting_activity = delivery_activity
    self.time_delivery = time_delivery
    self.applications = []

  def current_activity(self, rate_constant):
    """Returns the current activity based on the time and applied quantities"""
    current_activity = self.starting_activity
    start_time = self.time_delivery

    #this loop calculates the decrease in activities for applications and times between applications
    for application in self.applications:
      if application[0]<datetime.now():
        time_difference = application[0] - start_time
        time_difference_minutes = time_difference.total_seconds() / 60
        current_activity = current_activity*exp(-rate_constant*(time_difference_minutes))
        current_activity -= application[1]
        start_time = application[0]

    #Now account for decrease of activity from the last application
    current_time = datetime.now()
    if current_time < self.time_delivery:
      print("Current activity of this stock cannot be requested as it has not been delivered yet.")
      return 0
    time_difference = current_time - start_time
    time_difference_minutes = time_difference.total_seconds() / 60
    current_activity = current_activity*exp(-rate_constant*(time_difference_minutes))

    return current_activity

  def timepoint_activity(self, timepoint, rate_constant):
    """Returns the activity at certain timepoint based on the time and applied quantities"""
    #if we request the timepoint_activity befofe its time of the delivery the function returns nothing
    if timepoint < self.time_delivery:
      #print("The timepoint activity can be requested only for the times after the stock has been delivered")
      return 0
    current_activity = self.starting_activity
    start_time = self.time_delivery
    #if there is no application yet at the desired timepoint loop cannot run
    if len(self.applications) == 0:
      current_time = timepoint
      time_difference = current_time - start_time
      time_difference_minutes = time_difference.total_seconds() / 60
      current_activity = current_activity * exp(-rate_constant * (time_difference_minutes))

    for application in self.applications:
      if application[0] < timepoint:
        time_difference = application[0] - start_time
        time_difference_minutes = time_difference.total_seconds() / 60
        current_activity = current_activity*exp(-rate_constant*(time_difference_minutes))
        current_activity -= application[1]
        start_time = application[0]
      else:
        time_difference = timepoint - start_time
        time_difference_minutes = time_difference.total_seconds() / 60
        current_activity = current_activity*exp(-rate_constant*(time_difference_minutes))

    return(current_activity)

#stores data for multiple different deliveries of 1 radiofarmacum
#deliveries variable in init is list of tuples containing the delivery time and delivery activity
class RadiofarmacumStock:
  def __init__(self, halftime, type_, deliveries):
    """Halftime input in minutes, starting activity in MBq"""
    self.halftime = halftime
    self.rate_constant = log(2)/halftime
    self.type_ = type_
    self.deliveries = [Delivery(time_delivery, delivery_activity) for time_delivery, delivery_activity in deliveries]

  def print_status(self):
    """This function prints out the current state of the stock"""
    deliveries_in_stock = [delivery for delivery in self.deliveries if delivery.time_delivery<datetime.now()]
    string = ""
    for delivery in deliveries_in_stock:
      string += f"The current activity of stock {self.type_} delivered at {delivery.time_delivery} is {delivery.current_activity(self.rate_constant)}\n"

    print(string)

  def application(self, application_time, applied_activity, id, delivery_time):
    """Takes the information of time and the activity of the applied treatment and stores it to the internal database of applications"""

    #selects the correct delivery, checks if the application time is after the delivery time and whether we have suffictient activity
    for delivery in self.deliveries:
      if delivery.time_delivery == delivery_time:
        if delivery.time_delivery <= application_time:
          if delivery.timepoint_activity(application_time, self.rate_constant) >= applied_activity:
            delivery.applications.append((application_time, applied_activity, id))
          else:
            print(f"There is not enough of activity for the radiofarmacum {self.type_}")
        else:
           print("The application cannot be done using this stock as it has not been delivered yet")

  def change_application_time(self, id, new_time, delivery_time):
    delivery_time_ = time_formatting(delivery_time)
    for delivery in self.deliveries:
      if delivery_time_ == delivery.time_delivery:
        for i, application in enumerate(delivery.applications):
          if application[2] == id:
            activity = application[1]
            delivery.applications.pop(i)
            self.application(new_time, activity, id, delivery_time)

  def timepoint_activity(self, timepoint, delivery_time):
    """Returns the activity at certain timepoint based on the time and applied quantities"""
    #if we request the timepoint_activity befofe its time of the delivery the function returns nothing
    delivery_time = time_formatting(delivery_time)
    timepoint = time_formatting(timepoint)
    for delivery in self.deliveries:
      if delivery.time_delivery == delivery_time:
        return delivery.timepoint_activity(timepoint, self.rate_constant)

  def current_activity(self, delivery_time):
    """Returns the current activity based on the time and applied quantities"""
    delivery_time = time_formatting(delivery_time)
    for delivery in self.deliveries:
      if delivery.time_delivery == delivery_time:
        return delivery.current_activity(self.rate_constant)

  def sum_current_activities(self):
    """Returns the current sum of activities of all stock samples"""
    sum_of_activities = 0
    for delivery in self.deliveries:
      sum_of_activities += delivery.current_activity(self.rate_constant)
    return sum_of_activities

  def sum_timepoint_activities(self, timepoint):
    """Returns the sum of activities of all stock samples at specified timepoint"""
    sum_of_activities = 0
    for delivery in self.deliveries:
      sum_of_activities += delivery.timepoint_activity(timepoint, self.rate_constant)
    return sum_of_activities

class Schedule:
    def __init__(self):
        self.treatments = {}  # patient_name -> treatment info

    def add_treatment(self, patient, application_time, accumulation_time_minutes, scan_duration_minutes):
        """Adds a treatment to the schedule."""
        application_time = self._format_minutes_datetime(application_time)
        accumulation = timedelta(minutes=accumulation_time_minutes)
        scan_duration = timedelta(minutes=scan_duration_minutes)
        scan_start = application_time + accumulation
        scan_end = scan_start + scan_duration
        self.treatments[patient] = {
            "application_time": application_time,
            "accumulation": accumulation,
            "scan_duration": scan_duration,
            "scan_start": scan_start,
            "scan_end": scan_end
        }

    def change_application_time(self, patient, new_application_time):
        """Updates the application time for the patient and recalculates PET scan times."""
        if patient not in self.treatments:
            print(f"No treatment found for patient {patient}")
            return
        new_time = self._format_time(new_application_time)
        self.treatments[patient]["application_time"] = new_time
        self._recalculate_scan_window(patient)

    def _recalculate_scan_window(self, patient):
        treatment = self.treatments[patient]
        treatment["scan_start"] = treatment["application_time"] + treatment["accumulation"]
        treatment["scan_end"] = treatment["scan_start"] + treatment["scan_duration"]

    def check_and_resolve_overlaps(self):
        """Detects and resolves scan time overlaps by postponing the minimum number of treatments."""
        # Sort by scan start
        sorted_patients = sorted(self.treatments.items(), key=lambda x: x[1]["scan_start"])
        for i in range(len(sorted_patients) - 1):
            current_patient, current = sorted_patients[i]
            next_patient, next_ = sorted_patients[i + 1]
            if current["scan_end"] > next_["scan_start"]:
                # Conflict! Move next patient's application time forward
                delay = current["scan_end"] - next_["scan_start"]
                new_application = next_["application_time"] + delay
                self.change_application_time(next_patient, new_application)
                # Re-sort and restart to make sure we preserve minimum changes
                return self.check_and_resolve_overlaps()

    def print_schedule(self):
        for patient, data in sorted(self.treatments.items(), key=lambda x: x[1]["scan_start"]):
            print(f"{patient}: Application at {data['application_time'].time()}, "
                  f"Scan from {data['scan_start'].time()} to {data['scan_end'].time()}")

    def _format_time(self, time_tuple):
        """Converts a time tuple (YYYY, MM, DD, hh, mm) to datetime object."""
        return datetime(*time_tuple)

    def _format_minutes_datetime(self, minutes):
       start_time = datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
       return start_time + timedelta(minutes=minutes)




#---- David fully funcitoning code - Daniel approved


class Person:
  def __init__(self, internal_idx, mass):
    self.internal_idx = internal_idx
    self.mass = mass




class RadiofarmacumType:
  def __init__(self, name, halftime, price_per_GBq, GBq_per_kg_or_dose, time_on_pet, time_accumulation, times_deliveries):
    self.name = name
    self.halftime = halftime
    self.rate_constant = math.log(2)/halftime
    self.GBq_per_kg_or_dose = GBq_per_kg_or_dose
    self.price_per_kg_or_dose = price_per_GBq * GBq_per_kg_or_dose
    self.time_on_pet = time_on_pet
    self.time_accumulation = time_accumulation
    self.times_deliveries = times_deliveries

  def calculate_base_price(self, mass):
    if (self.name == "vizamyl" or self.name == "FDG_brain"):
      return self.price_per_kg_or_dose
    return mass * self.price_per_kg_or_dose

  def calculate_price(self, mass, time_from_delivery):
    return self.calculate_base_price(mass) * math.exp(time_from_delivery * self.rate_constant)

  def set_time_delivery(self, time):
    self.time_delivery = time

  def get_time_delivery(self):
    return self.time_delivery

  def get_activity_needed_treatment(self, mass):
    if (self.name == "vizamyl" or self.name == "FDG_brain"):
      return self.GBq_per_kg_or_dose
    return mass * self.GBq_per_kg_or_dose

  def get_activity_needed_order(self, mass, time_acc):
    return self.get_activity_needed_treatment(mass) * math.exp((time_acc - self.get_time_delivery()) * self.rate_constant)


  def get_time_from_delivery(self, time):
    """Returns the time from last delivery in minutes, time in minutes from midnight"""
    smaller = [delivery for delivery in self.times_deliveries if delivery <= time]
    if len(smaller) == 0: return None
    return (time - max(smaller))




clinic_opening_time =  6 * 60 + 30 #(6h + 30mins)
FDG_onco  = RadiofarmacumType("FDG_onco",   halftime=109.8,  price_per_GBq = 8900,    GBq_per_kg_or_dose = 0.0025,  time_on_pet = 25, time_accumulation = 60,  times_deliveries = [6 * 60 + 30, 10 * 60 + 30, 13 * 60] )
FDG_brain = RadiofarmacumType("FDG_brain",  halftime=109.8,  price_per_GBq = 8900,    GBq_per_kg_or_dose = 0.15,    time_on_pet = 60, time_accumulation = 0,   times_deliveries = [6 * 60 + 30, 10 * 60 + 30, 13 * 60] )
Ga68      = RadiofarmacumType("Ga68",       halftime=68,     price_per_GBq = 10_000,  GBq_per_kg_or_dose = 0.002,   time_on_pet = 30, time_accumulation = 60,  times_deliveries = [7 * 60, 11 * 60])
C11       = RadiofarmacumType("methionin",  halftime=20.4,   price_per_GBq = 25_000,  GBq_per_kg_or_dose = 0.0045,  time_on_pet = 20, time_accumulation = 90,  times_deliveries = [10 * 60])
F18       = RadiofarmacumType("vizamyl",    halftime=109.8,  price_per_GBq = 75_000,  GBq_per_kg_or_dose = 0.185,   time_on_pet = 20, time_accumulation = 90,  times_deliveries = [11 * 60])
pharms_dict = {"FDG_onco": FDG_onco, "FDG_brain": FDG_brain, "Ga68": Ga68, "methionin": C11, "vizamyl": F18}



def generate_test_set(size):
  return [(Person(i, i, rnd.gauss(85.0, 20.0)), pharms_dict[rnd.randint(0, len(pharms_dict) - 1)]) for i in range(size)]


def greatest_less_than(lst, x):
    smaller = [num for num in lst if num <= x]
    return max(smaller)

def calculate_starting_times(patient_order):
  """Return (pet_times, inj_times)"""
  time_pet = []
  time_inj = []
  tot_time = clinic_opening_time
  for (person, pharm) in patient_order:
    if (tot_time < pharm.time_accumulation + clinic_opening_time):
      tot_time = clinic_opening_time + pharm.time_accumulation
    if (pharm.get_time_from_delivery(tot_time) is None):
      tot_time = pharm.times_deliveries[0] + pharm.time_accumulation
    time_pet.append(tot_time)
    time_inj.append(tot_time - pharm.time_accumulation)
    tot_time += pharm.time_on_pet
  return (time_pet, time_inj)

def evaluate_solution(order, pet_time_starts):
  cost = 0
  for i in range(len(order)):
      (person, pharm) = order[i]
      inj_time_start = pet_time_starts[i] - pharm.time_accumulation
      time_from_delivery_to_injection = pharm.get_time_from_delivery(inj_time_start)
      if time_from_delivery_to_injection is None: return None
      cost += pharm.calculate_price(person.mass, time_from_delivery_to_injection)
      # print(f"Cost: {pharm.calculate_price(person.mass, time_from_delivery_to_pet)}")
      # print(f"Person {person.id}, mass {person.mass}")
      # print(f"Radiopharm {pharm.name}, halftime {pharm.rate_constant}, $/kg {pharm.price_per_kg_or_dose}, time {pharm.time_on_pet}")
      # print(f"Time starts: {time_starts[i]}")
  return cost

def find_optimum(test_set):
  lst = [i for i in range(len(test_set))]
  perms = [list(p) for p in itertools.permutations(lst)]
  best_cost = math.inf
  best_order = None
  best_time_pet = None
  best_time_inj = None
  for p in perms:
    # print("---------------")
    patient_order = [test_set[i] for i in p]
    (times_starts_pet, times_starts_inj) = calculate_starting_times(patient_order)
    if times_starts_pet is None: continue
    cost = evaluate_solution(patient_order, times_starts_pet)
    if cost is None: continue
    # print(f"Time pet starts: {times_starts_pet}")
    # print(f"Time inj starts: {times_starts_inj}")
    # print(f"Time deliveries: {time_deliveries}")
    # print(f"Total cost: {cost}")
    if cost < best_cost:
      best_cost = cost
      best_order = patient_order
      best_time_pet = times_starts_pet
      best_time_inj = times_starts_inj
  return (best_cost, best_order, best_time_pet, best_time_inj)

# (best_cost, best_order, best_time_pet, best_time_inj) = find_optimum(test_set)
# print(f"Best price: {best_cost}")
# for (person, pharm) in best_order:
#   print(f"Person {person.id}, mass {person.mass}")
#   print(f"Radiopharm {pharm.name}, halftime {pharm.halftime}, $/kg {pharm.price_per_kg_or_dose}, time {pharm.time_accumulation}")
#   # time_from_delivery = pharm.get_time_from_delivery()
#   # print(f"Cost of the person {pharm.calculate_price()}")
# print(f"Start pet {best_time_pet}")
# print(f"Start inj {best_time_inj}")




def constraint_programming(test_set):
  N = len(test_set)
  SLOT_DURATION = 5  # minut
  delivery_slots = [
      [int(math.ceil((t - clinic_opening_time)/ SLOT_DURATION))
       for t in pharm.times_deliveries] for (_, pharm) in test_set
  ]
  duration = [int(math.ceil(pharm.time_on_pet / SLOT_DURATION)) for (_, pharm) in test_set]
  application_slots = [int(math.ceil(pharm.time_accumulation / SLOT_DURATION)) for (_, pharm) in test_set]
  base_prices = [int(pharm.calculate_base_price(person.mass)) for (person, pharm) in test_set]
  decay_rates = [pharm.rate_constant for (_, pharm) in test_set]  # rychlost rozpadu

  # Model
  model = cp_model.CpModel()

  # Max. slot
  max_slot = 120  # = 12 hodin * 60 / 5

  # Start time proměnné
  start_times = [model.NewIntVar(0, max_slot - duration[i], f"start_{i}") for i in range(N)]

  # Žádné překryvy
  for i in range(N):
      for j in range(i + 1, N):
          is_earlier_i = model.NewBoolVar(f"is_earlier_{i}_{j}")
          is_earlier_j = model.NewBoolVar(f"is_earlier_{j}_{i}")

          model.Add(start_times[i] + duration[i] <= start_times[j]).OnlyEnforceIf(is_earlier_i)
          model.Add(start_times[j] + duration[j] <= start_times[i]).OnlyEnforceIf(is_earlier_j)

          model.AddBoolOr([is_earlier_i, is_earlier_j])

  # Celkové náklady
  total_cost = []
  for i in range(N):
      # Najdi proměnnou: nejbližší dovoz <= start
      delay_costs = []

      for t in range(max_slot + 1):
          # Najdi nejbližší dovoz
          if t < application_slots[i]:
            delay_costs.append(10_000_000)  # or some max cost / penalty
            continue

          try:
            latest_delivery = max(d for d in delivery_slots[i] if d <= t - application_slots[i])
          except ValueError:
            delay_costs.append(10_000_000)
            continue

          delay_minutes = (t - application_slots[i] - latest_delivery) * SLOT_DURATION
          price = int(math.exp(decay_rates[i] * delay_minutes) * base_prices[i])
          delay_costs.append(price)

      max_cost = max(delay_costs) if delay_costs else 10_000_000
      cost_var = model.NewIntVar(0, max_cost, f"cost_{i}")
      model.AddElement(start_times[i], delay_costs, cost_var)
      total_cost.append(cost_var)

  # Minimalizace
  model.Minimize(sum(total_cost))

  # Solver
  solver = cp_model.CpSolver()
  status = solver.Solve(model)

  # Výpis výsledků
  if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
      # print('Řešení nalezeno!')
      lst = []
      for i in range(N):
          start = solver.Value(start_times[i]) * SLOT_DURATION + clinic_opening_time
          lst.append((i, start))
          # print(f"Pacient {i}: start v {start * Delka_slotu} min, konec v {(start + duration[i]) * Delka_slotu} min")
      sorted_lst = sorted(lst, key=lambda x: x[1])
      results = [(test_set[i], start) for (i,start) in sorted_lst]
      print(f"Celkové náklady: {solver.ObjectiveValue():.2f}")
      return results
  else:
      print("Žádné řešení.")


class ResultData:
    def __init__(self, best_order_and_times):
        for ((person, pharm), time_acc) in best_order_and_times:
            time_of_delivery = time_acc - pharm.time_accumulation
            latest_delivery = max(d for d in pharm.times_deliveries if d <= time_of_delivery)
            pharm.set_time_delivery(latest_delivery)
        self.best_order_and_times = sorted(best_order_and_times, key=lambda x: x[0][0].internal_idx)
        self.itinerary = self.create_itinerary()

    def create_itinerary(self):
        itinerary_data = []
        for ((person, pharm), time) in self.best_order_and_times:
            # event, who, time (mins from midnight)
            one_event_pet = ("pet_starts", person.internal_idx, time)
            itinerary_data.append(one_event_pet)
            one_event_acc = ("acc_starts", person.internal_idx, time - pharm.time_accumulation)
            itinerary_data.append(one_event_acc)
        for name, radiofarmacum in pharms_dict.items():
            for time in radiofarmacum.times_deliveries:
                itinerary_data.append(("new_pharm", name, time))

        sorted_itinerary_data = sorted(itinerary_data, key=lambda x: x[2])
        return sorted_itinerary_data

    def sort_by_lowest_accumulation_time(self):
        self.best_order_and_times = sorted(self.best_order_and_times, key=lambda x: x[1])
    def sort_by_idx(self):
        self.best_order_and_times = sorted(self.best_order_and_times, key=lambda x: x[0][0].internal_idx)
    def get_starting_accumulation_times(self):
        return [time - pharm.time_accumulation for ((person, pharm), time) in self.best_order_and_times]
    def get_total_pet_times(self):
        return [pharm.time_on_pet for ((person, pharm), time) in self.best_order_and_times]
    def get_total_accumulation_times(self):
        return [pharm.time_accumulation for ((person, pharm), time) in self.best_order_and_times]
    def get_constant_rate(self):
        return [pharm.rate_constant() for ((person, pharm), time) in self.best_order_and_times]
    def get_activity_needed_treatment(self):
        return [pharm.get_activity_needed_treatment(person.mass) for ((person, pharm), time) in self.best_order_and_times]
    def get_activity_needed_order(self):
        return [pharm.get_activity_needed_order(person.mass, time - pharm.time_accumulation) for ((person, pharm), time) in self.best_order_and_times]
    def get_order_data(self):
        order_data = []
        for ((person, pharm), time) in self.best_order_and_times:
            one_order = {
                "type" : pharm.name,
                "id" : person.internal_idx,
                "delivery_time" : pharm.get_time_delivery(),
                "activity_order" : pharm.get_activity_needed_order(person.mass, time - pharm.time_accumulation),
                "activity_treatment" : pharm.get_activity_needed_treatment(person.mass),
                "application_time" : time - pharm.time_accumulation,
                "half_time" : pharm.halftime
            }
            order_data.append(one_order)
        return order_data
    def get_itinerary(self, daytime):
        """ Returns the list of 2 tuples (what, who, time_from_now)"""
        minutes = daytime.hour * 60 + daytime.minute
        filtered_itinerary = [item for item in self.itinerary if item[2] >= minutes][:2]
        midnight = datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
        return [(item[0], item[1], item[2] - minutes) for item in filtered_itinerary]

def generate_results(): #patients_data[idx] = {"internal_id":idx+1, "mass": mass, "id": patient_id, "treatment": treatment, "diabetes": diabetes}
    test_set = generate_test_set(10)
    best_order_and_times = constraint_programming(test_set)
    results = ResultData(best_order_and_times)
    return results


def generate_results_param(parameters):
    test_set = []
    for param in parameters:
        person = Person(internal_idx=param["id_entry"], mass=param["mass_entry"])
        pharm = copy.deepcopy(pharms_dict[param["treatment_var"]])
        test_set.append([person, pharm])
    # test_set = generate_test_set(10)
    best_order_and_times = constraint_programming(test_set)
    results = ResultData(best_order_and_times)
    return results



#---- Displaying results


# Global variables
number_of_patients = 0
patients_data = []
start_times = []
radiofarmacum_stock = []
PET_times = []
accumulation_times = []

def Generate_patients(optimil_func_output):
   schedule = Schedule()
   for pat in optimil_func_output:
      schedule.add_treatment(pat["id"], pat["application_time"], pat["accumulation_time"], pat["scan_time"])


# Function to open a new window with the form for patient application and dose entry
def open_dose_form():
    # Create a new window
    dose_window = tk.Toplevel(root)
    dose_window.title("Patient Dose Application")
    dose_window.geometry("600x400")  # Adjust the window size as needed

    tk.Label(dose_window, text="Patient Dose Application Form", font=("Helvetica", 16)).pack(pady=10)

    # Frame for checkboxes and dose fields
    form_frame = tk.Frame(dose_window)
    form_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    # Create checkboxes and entry fields for each patient
    applied_doses = []  # List to store dose data
    checkboxes = []  # List to store checkbox variables for patient selection

    for idx, patient in enumerate(patients_data):
        frame = tk.Frame(form_frame)
        frame.pack(anchor="w", pady=5)

        # Patient selection checkbox
        var = tk.BooleanVar(value=False)  # Boolean variable for the checkbox
        tk.Checkbutton(frame, text=f"Patient {patient['id_entry']} (Mass: {patient['mass_entry']}kg)", variable=var, font=("Helvetica", 12)).pack(side=tk.LEFT)
        checkboxes.append(var)

        # Dose entry field
        tk.Label(frame, text="Applied Dose (mBq):", font=("Helvetica", 12)).pack(side=tk.LEFT, padx=5)
        dose_entry = tk.Entry(frame, font=("Helvetica", 12))
        dose_entry.pack(side=tk.LEFT)
        applied_doses.append(dose_entry)

    # Submit button
    def handle_submit():
        # Collect data from the form
        results = []
        for idx, checkbox in enumerate(checkboxes):
            if checkbox.get():  # If the patient is selected (checkbox is checked)
                dose = applied_doses[idx].get()
                if not dose:
                    tk.messagebox.showerror("Error", f"Please enter a dose for Patient {patients_data[idx]['id']}.")
                    return
                try:
                    dose = float(dose)  # Validate dose as a numerical value
                except ValueError:
                    tk.messagebox.showerror("Error", f"Invalid dose value for Patient {patients_data[idx]['id']}.")
                    return
                results.append({
                    "Patient ID": patients_data[idx]["id"],
                    "Dose Applied (mBq)": dose
                })

        if not results:
            tk.messagebox.showinfo("No Selection", "No patients were selected for application.")
        else:
            # Display results (you can replace this with another action, e.g., saving data)
            print("Applied Doses:")
            for entry in results:
                print(entry)
            tk.messagebox.showinfo("Success", "Data submitted successfully!")

    submit_button = tk.Button(dose_window, text="Submit", font=("Helvetica", 12), bg="green", fg="white", command=handle_submit)
    submit_button.pack(pady=20)




# Function for the button action
def open_activity_window():
    # Create a new window
    activity_window = tk.Toplevel(root)
    activity_window.title("Radiopharmaceuticals Activity")
    #activity_window.geometry("1000x600")  # Set window size

    tk.Label(activity_window, text="Grafy aktivit radiofarmak", font=("Helvetica", 16)).pack(pady=10)

    # Frame to hold the graphs
    graph_frame = tk.Frame(activity_window)
    graph_frame.pack(fill=tk.BOTH, expand=True)





    # Number of unique treatments
    num_radiofarm = len(radiofarmacum_stock)

        # Create subplots for stacked graphs sharing the same x-axis

    fig, axs = plt.subplots(
        num_radiofarm, 1, figsize=(8, 6), sharex=True, gridspec_kw={'hspace': 0}
        )  # No space between graphs

        # Ensure axs is iterable, even if there's only one treatment
    if num_radiofarm == 1:
        axs = [axs]

    start_day = datetime.today().replace(hour=6, minute=0, second=0, microsecond=0)
    end_day = datetime.today().replace(hour=15, minute=30, second=0, microsecond=0)
    x = [start_day + timedelta(seconds=5*i) for i in range(round(((12*60+30)*60)))]

        # Generate data and plot for each treatment
    for idx, object in enumerate(radiofarmacum_stock):
        activities = []
        for i in range(0, len(x)):
            activities.append(object.sum_timepoint_activities(x[i]))


        axs[idx].plot(x, activities, label=object.type_, color=f"C{idx}")
            #axs[idx].set_title(f"{treatment} Aktivita", fontsize=40)


            # Configure x-axis format to match the main graph
        now = datetime.now()
        axs[idx].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        axs[idx].xaxis.set_major_locator(mdates.HourLocator(interval=1))  # Show hours
        axs[idx].tick_params(axis="x", which="both", bottom=False, labelbottom=(idx == num_radiofarm - 1))  # Bottom labels only for the last graph
        axs[idx].set_xlim(now - timedelta(minutes=60), now + timedelta(minutes=240))
        axs[idx].legend(object.type_)
        print(object.type_)
        custom_ticks = np.linspace(min(activities), max(activities)+np.average(activities), 4)
        axs[idx].set_yticks(custom_ticks)

        if idx == num_radiofarm - 1:  # Add xlabel to the last graph
            axs[idx].set_xlabel("Čas", fontsize=20)

    def update_red_line():
        now = datetime.now()
        ax.axvline(now, color='red', linestyle='--', linewidth=2, alpha=0.7)  # Red line for current time
        ax.axvspan(now - timedelta(minutes=5), now + timedelta(minutes=5), color='gray', alpha=0.05)
        axs[0].text(now + timedelta(minutes=10), 5,now.strftime("%H:%M"),color='red', ha='left', va='bottom', fontsize= 15, fontweight='bold')

        fig.canvas.draw()
        activity_window.after(1000, update_red_line)  # Redraw every second

    update_red_line()  # Initialize the red line update
        # Embed the graph into the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=graph_frame)
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    canvas.draw()



# Function to compute intervals
def get_intervals():
    today = datetime.today().replace(hour=0, minute=0, second=0, microsecond=0)
    intervals = []

    for idx, patient in enumerate(patients_data):
        start = today + timedelta(minutes=start_times[idx])
        accumulation_end = start + timedelta(minutes=accumulation_times[idx])
        pet_end = accumulation_end + timedelta(minutes=PET_times[idx])
        intervals.append((
            "Pacient" + '\n' + f"R.č.: {patient['id_entry']}",
            start,
            accumulation_end,
            pet_end
        ))
    return intervals

# Function to dynamically generate patient fields when number of patients is updated
def update_patient_fields():
    for widget in patient_fields_frame.winfo_children():
        widget.destroy()
    patients_data.clear()

    try:
        global number_of_patients
        number_of_patients = int(num_patients_entry.get())
        if number_of_patients <= 0:
            raise ValueError("Number of patients must be positive.")

        for i in range(number_of_patients):
            frame = tk.Frame(patient_fields_frame)
            frame.pack(pady=5, anchor="w")

            tk.Label(frame, text=f"Pacient {i + 1} Hmotnost (kg):", font=("Helvetica", 12)).grid(row=0, column=0)
            mass_entry = tk.Entry(frame, font=("Helvetica", 12))
            mass_entry.grid(row=0, column=1)

            tk.Label(frame, text="Rodné číslo:", font=("Helvetica", 12)).grid(row=0, column=2)
            id_entry = tk.Entry(frame, font=("Helvetica", 12))
            id_entry.grid(row=0, column=3)

            tk.Label(frame, text="Typ léčby:", font=("Helvetica", 12)).grid(row=0, column=4)
            treatment_var = tk.StringVar(value="FDG_onco")
            treatment_menu = tk.OptionMenu(frame, treatment_var, "FDG_onco", "FDG_brain", "Ga68", "methionin", "vizamyl")
            treatment_menu.grid(row=0, column=5)

            diabetes_var = tk.BooleanVar(value=False)
            diabetes_checkbox = tk.Checkbutton(frame, text="Diabetes", variable=diabetes_var, font=("Helvetica", 12))
            diabetes_checkbox.grid(row=0, column=6)

            patients_data.append({"mass_entry": mass_entry, "id_entry": id_entry, "treatment_var": treatment_var, "diabetes_var": diabetes_var})

    except ValueError:
        tk.messagebox.showerror("Invalid Input", "Please enter a valid number of patients.")

# Function to handle form submission
def handle_form_submission():
    global start_times
    global PET_times
    global accumulation_times
    global radiofarmacum_stock
    global results

    try:
        for idx, patient in enumerate(patients_data):
            mass = float(patient["mass_entry"].get())
            patient_id = patient["id_entry"].get()
            treatment = patient["treatment_var"].get()
            diabetes = patient["diabetes_var"].get()

            if not patient_id.strip():
                raise ValueError("ID cannot be empty.")

            patients_data[idx] = {"mass_entry": mass, "id_entry": patient_id, "treatment_var": treatment, "diabetes_var": diabetes}


        results = generate_results_param(patients_data)
        results.sort_by_lowest_accumulation_time()
        start_times = results.get_starting_accumulation_times()    #[480]*number_of_patients   #? Here function from David for getting start times, input: patients_data, output: dict{id, start_times}
        PET_times = results.get_total_pet_times()
        accumulation_times = results.get_total_accumulation_times()

        list_radiofarmacum_patient_needs = results.get_order_data()
        dictionary_of_radiofarmacums = {}
        #time today midnight
        midnight = time_formatting((2025, 4, 13, 0, 0))
        for patient in list_radiofarmacum_patient_needs:
            type = patient["type"]
            activity_order = patient["activity_order"]
            delivery_time = midnight + timedelta(minutes=patient["delivery_time"])
            halftime = patient["half_time"]
            if type not in dictionary_of_radiofarmacums:
                dictionary_of_radiofarmacums[type] = (halftime, {})
                if delivery_time not in dictionary_of_radiofarmacums[type][1]:
                    dictionary_of_radiofarmacums[type][1][delivery_time] = activity_order
                else:
                    dictionary_of_radiofarmacums[type][1][delivery_time] += activity_order
            else:
                if delivery_time not in dictionary_of_radiofarmacums[type]:
                    dictionary_of_radiofarmacums[type][1][delivery_time] = activity_order
                else:
                    dictionary_of_radiofarmacums[type][1][delivery_time] += activity_order

        for radiofarmacum in dictionary_of_radiofarmacums:
            deliveries = []
            radiofarmacum_halftime = dictionary_of_radiofarmacums[radiofarmacum][0]
            for delivery in dictionary_of_radiofarmacums[radiofarmacum][1]:
                deliveries.append((delivery, dictionary_of_radiofarmacums[radiofarmacum][1][delivery]))
            radiofarmacum_stock.append(RadiofarmacumStock(radiofarmacum_halftime, radiofarmacum, deliveries))

        for patient in list_radiofarmacum_patient_needs:
            time_application = midnight + timedelta(minutes=patient["application_time"])
            time_delivery = midnight + timedelta(minutes=patient["delivery_time"])
            id = patient["id"]
            activity_treatment = patient["activity_treatment"]
            for farmac in radiofarmacum_stock:
                if patient["type"] == farmac.type_:
                    farmac.application(time_application, activity_treatment, id, time_delivery)

        #test not really needed
        # for rad in radiofarmacum_stock:
        #     print(rad.type_, rad.halftime, len(rad.deliveries),[(delivery.time_delivery, delivery.starting_activity) for delivery in rad.deliveries])

        #Generate_patients(patients_data)

        form_frame.pack_forget()
        display_graph()

    except ValueError as e:
        tk.messagebox.showerror("Invalid Input", str(e))

# Function to display the graph
def display_graph():
    graph_frame.pack(fill=tk.BOTH, expand=True)
    update_chart()

# Function to update the graph dynamically
def update_chart():
    draw_chart(ax)
    canvas.draw()
    root.after(1000, update_chart)  # Redraw every second to update the red time line

# Function to draw the chart
def draw_chart(ax):
    ax.clear()
    intervals = get_intervals()

    start_day = datetime.today().replace(hour=6, minute=0, second=0, microsecond=0)
    end_day = datetime.today().replace(hour=19, minute=0, second=0, microsecond=0)

    yticks = []
    ylabels = []

    for i, (label, start, acc_end, pet_end) in enumerate(intervals):
        y_pos = i * 2
        if i%2 == 1:
            color = "blue"
        else:
            color = "#ADF2A7"

        ax.barh(y_pos, acc_end - start, left=start, height=2, color=color, alpha=0.5)  # Accumulation
        ax.barh(y_pos, pet_end - acc_end, left=acc_end, height=2, color=color)  # PET scan

        yticks.append(y_pos)
        ylabels.append(label)

    ax.set_xlim(start_day, end_day)
    ax.set_ylim(-1, number_of_patients * 2)
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=15 - round(number_of_patients/1.5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_title("Rozvrh PET vyšetření", fontsize=25)

    now = datetime.now()
    ax.axvline(now, color='red', linestyle='--', linewidth=5)  # Red time now line
    ax.axvspan(now - timedelta(minutes=5), now + timedelta(minutes=5), color='gray', alpha=0.2)
    ax.text(
        now + timedelta(minutes=10),
        number_of_patients * 2 - 2*(number_of_patients/7),
        now.strftime("%H:%M"),
        color='red', ha='left', va='bottom', fontsize=20 - round(number_of_patients/4), fontweight='bold'
    )

    minor_ticks = []
    for i in range(number_of_patients):
        y_pos = i * 2
        bar_height = 2
        minor_ticks.extend([y_pos - bar_height / 2, y_pos + bar_height / 2])

    ax.set_yticks(minor_ticks, minor=True)
    ax.tick_params(axis='y', which='minor', length=5, color='gray')
    ax.yaxis.grid(True, which='minor', linestyle='--', linewidth=0.7, alpha=1)


    # def make_text():
    #     next_actions = results.get_itinerary(now)
    #     top_right_text = ax.text(ax.get_xlim()[1]-0.02, ax.get_ylim()[1]-0.02,  "", ha="right", va="top", fontsize = 16, color="blue")

    #     text = ""
    #     for (what, who, time) in next_actions:
    #         if str(what) == "pet_starts":
    #             matter = "Scanning začíná"  + " pacient " + str(who)
    #         elif str(what) == "acc_starts":
    #            matter = "Akumulaci začíná"  + " pacient " + str(who)
    #         else:
    #             if str(who) == "FDG_onco" or "FDG_brain":
    #               matter = "Dodávka " + "FDG"
    #             else:
    #                matter = "Dodávka " + str(who)

    #         text += str(time) + " min: \n" +  str(matter) + "\n"

    #     top_right_text.set_text(text)
    def make_text():
      next_actions = results.get_itinerary(now)

      # Initialize an empty string for the text
      text = ""
      for (what, who, time) in next_actions:
          if str(what) == "pet_starts":
              matter = "Scanning začíná" + " pacient " + str(who)
          elif str(what) == "acc_starts":
              matter = "Akumulaci začíná" + " pacient " + str(who)
          else:
              if str(who) == "FDG_onco" or str(who) == "FDG_brain":
                  matter = "Dodávka " + "FDG"
              else:
                  matter = "Dodávka " + str(who)

          text += "\n" + str(time) + " min: \n" + str(matter)

      # Add the text in the top-right corner (placeholder position)
      top_right_text = ax.text(
          ax.get_xlim()[1] - 0.02, ax.get_ylim()[1] - 0.03,
          text, ha="right", va="top", fontsize=16, color="blue"
      )

      # Get the bounding box of the text
      bbox = top_right_text.get_window_extent(renderer=ax.figure.canvas.get_renderer())
      bbox_data_coords = ax.transData.inverted().transform(bbox)  # Transform bbox to data coordinates

      # Extract coordinates for the rectangle
      x_min, y_min = bbox_data_coords[0]  # Bottom-left corner
      x_max, y_max = bbox_data_coords[1]  # Top-right corner

      # Extend the box in the x-direction for padding
      padding = 0.1 * (x_max - x_min)  # Add 10% extra width as padding
      x_min -= padding
      x_max += padding

      # Draw a dashed rectangle around the text
      rect = patches.Rectangle(
          (x_min, y_min),  # Bottom-left corner
          x_max - x_min,  # Width
          y_max - y_min,  # Height
          linewidth=1,
          edgecolor="red",
          linestyle="--",
          facecolor="none"
      )
      ax.add_patch(rect)  # Add rectangle to the axis

    make_text()
    #plt.tight_layout()

# GUI setup
root = tk.Tk()
root.title("PET Treatment Scheduler")

form_frame = tk.Frame(root)
form_frame.pack(fill=tk.BOTH, expand=True)

tk.Label(form_frame, text="Počet pacientů na den", font=("Helvetica", 14)).pack(pady=10)
num_patients_entry = tk.Entry(form_frame, font=("Helvetica", 12))
num_patients_entry.pack(pady=5)
num_patients_entry.bind("<Return>", lambda event: update_patient_fields())

patient_fields_frame = tk.Frame(form_frame)
patient_fields_frame.pack(pady=10)

generate_fields_button = tk.Button(form_frame, text="Generuj políčka pro pacienty", font=("Helvetica", 12),
                                    command=update_patient_fields, bg="blue", fg="white")
generate_fields_button.pack(pady=10)

submit_button = tk.Button(form_frame, text="Potvrď", font=("Helvetica", 12),
                           command=handle_form_submission, bg="green", fg="white")
submit_button.pack(pady=20)


graph_frame = tk.Frame(root)
fig, ax = plt.subplots(figsize=(16, 9))
canvas = FigureCanvasTkAgg(fig, master=graph_frame)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)


live_label = tk.Label(
    graph_frame,
    text="Živě",
    bg="red",
    fg="white",
    font=("Helvetica", 16),
    padx=10,
    pady=5
)
live_label.place(relx=0.95, rely=0, anchor="ne")  # Position on the top right corner

# View Activity button
view_button = tk.Button(
    graph_frame,
    text="Aktivita radiofarmak",
    font=("Helvetica", 16),
    bg="blue",
    fg="white",
    command=open_activity_window,
    width=16,  # Matches size
    height=1,  # Matches padding of the Live button
    relief="groove"
)
view_button.place(relx=0.9, rely=0.99, anchor="se")  # Positioned slightly left of the Live button



dose_button = tk.Button(
    graph_frame,
    text="Aplikace",
    font=("Helvetica", 16),
    bg="blue",
    fg="white",
    command=open_dose_form,
    width=16,
    height=1,
    relief="groove"
)
dose_button.place(relx=0.125, rely=0.99, anchor="sw")  # Position below the graph on the left



root.mainloop()