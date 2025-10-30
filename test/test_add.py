import unittest

def add(a, b):
    """Simple function to add two numbers."""
    return a + b


class TestAdd(unittest.TestCase):
    """Test cases for the add function."""
    
    def test_add_positive_numbers(self):
        """Test adding two positive numbers."""
        self.assertEqual(add(2, 3), 5)
        self.assertEqual(add(10, 20), 30)
    
    def test_add_negative_numbers(self):
        """Test adding two negative numbers."""
        self.assertEqual(add(-5, -3), -8)
        self.assertEqual(add(-10, -15), -25)
    
    def test_add_mixed_numbers(self):
        """Test adding positive and negative numbers."""
        self.assertEqual(add(5, -3), 2)
        self.assertEqual(add(-10, 15), 5)
    
    def test_add_with_zero(self):
        """Test adding with zero."""
        self.assertEqual(add(0, 5), 5)
        self.assertEqual(add(10, 0), 10)
        self.assertEqual(add(0, 0), 0)
    
    def test_add_floats(self):
        """Test adding floating point numbers."""
        self.assertAlmostEqual(add(2.5, 3.7), 6.2)
        self.assertAlmostEqual(add(1.1, 2.2), 3.3)


if __name__ == '__main__':
    unittest.main()