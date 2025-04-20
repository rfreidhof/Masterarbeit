#This is a small test written to verify that events work properly.

import copy
import simulate_asgs
import unittest



class TestEvent(unittest.TestCase):
    def test_event(self):
        theta = 0.5
        R = 10
        vs = [0.3, 0.7]
        sigmas = [3,0]
        for i in range(1000):

            L_R, eventlist = infinite.generate_graph(sigmas, theta, R, vs)
            ancestors = infinite.generate_potential_ancestors(L_R, theta, vs, sigmas, True)
            individuals = copy.deepcopy(ancestors)


            for event in reversed(eventlist):

                match event:
                    case infinite.MutationEvent():
                        print("Attempting mutation test")
                        prev_type = individuals[event.target].typ
                        event.apply_event_forward(individuals)

                        self.assertEqual(individuals[event.target].typ, event.typ, "Mutations verändert Typ nicht korrekt")
                        self.assertEqual(prev_type, event.previousType, "Mutations speichert previoust type nicht korrekt")

                    case infinite.CoalescingEvent():
                        print("Attempting coalescing test")
                        prev_length = len(individuals)
                        event.apply_event_forward(individuals)

                        self.assertEqual(len(individuals)-1, prev_length, "Coalescence does not create a new individual.")
                        self.assertEqual(individuals[event.origin].typ, individuals[event.target].typ, "Coalescence type is not properly inherited.")
                        self.assertEqual(individuals[event.origin].ancestor_type, individuals[event.target].ancestor_type, "Ancestor is not properly inherited.")

                    case infinite.BranchingEvent():
                        print("Attempting branching test")
                        prev_length = len(individuals)
                        prev_target = individuals[event.target]
                        origin = individuals[event.origin]
                        event.apply_event_forward(individuals)

                        self.assertEqual(len(individuals)+1, prev_length, "Branching does not remove an individual.")
                        self.assertEqual(len(individuals), event.origin, "Branching does not remove the last individual.")
                        if origin.typ <= event.typ:
                            self.assertEqual(origin.typ, individuals[event.target].typ, "Branching type is not properly inherited.")
                            self.assertEqual(origin.ancestor_type, individuals[event.target].ancestor_type, "Branching ancestor is not properly inherited.")
                        else:
                            self.assertEqual(individuals[event.target], prev_target, "Silent branching changes something.")






if __name__ == '__main__':
    unittest.main()