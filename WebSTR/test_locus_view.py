import unittest
from unittest.mock import patch, MagicMock
import requests
import pandas as pd
from locus_view import *
import logging

# Set up logging for debugging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

class TestLocusView(unittest.TestCase):

    @patch('locus_view.requests.get')
    def test_GetSTRMetadataAPI_success(self, mock_get):
        """Test GetSTRMetadataAPI when the API returns valid data."""
        logging.info("Starting test: test_GetSTRMetadataAPI_success")
        
        # Sample mock response data
        mock_response_data = {
            "chr": "chr4",
            "start": 123456,
            "end": 123789,
            "motif": "AAG",
            "copies": 5,
            "gene_name": "GENE1",
            "gene_desc": "Description",
            "total_calls": 50,
            "frac_variable": 0.5,
            "avg_size_diff": 2.3
        }
        logging.debug(f"Mock response data: {mock_response_data}")
        
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = mock_response_data

        # Call the function
        repeat_id = "123"
        result = GetSTRMetadataAPI(repeat_id)

        # Assertions and additional logging
        logging.debug(f"Result from GetSTRMetadataAPI: {result}")
        self.assertEqual(result["chrom"], "chr4")
        self.assertEqual(result["start"], 123456)
        self.assertEqual(result["end"], 123789)
        self.assertEqual(result["motif"], "AAG")
        self.assertEqual(result["copies"], 5)
        self.assertEqual(result["crc_data"], [50, 0.5, 2.3])
        
        logging.info("Completed test: test_GetSTRMetadataAPI_success")

    @patch('locus_view.requests.get')
    def test_GetSTRMetadataAPI_failure(self, mock_get):
        """Test GetSTRMetadataAPI when the API returns an error."""
        logging.info("Starting test: test_GetSTRMetadataAPI_failure")
        
        # Simulate an error response from the API
        mock_get.return_value.status_code = 404
        mock_get.return_value.json.return_value = {}

        repeat_id = "123"
        result = GetSTRMetadataAPI(repeat_id)
        
        # Since the API call failed, result should be None
        logging.debug(f"Result from GetSTRMetadataAPI on failure: {result}")
        self.assertIsNone(result)

        logging.info("Completed test: test_GetSTRMetadataAPI_failure")

    @patch('locus_view.requests.get')
    def test_GetFreqSTRInfoAPI_success(self, mock_get):
        """Test GetFreqSTRInfoAPI when API returns frequency data."""
        logging.info("Starting test: test_GetFreqSTRInfoAPI_success")
        
        # Sample mock response data
        mock_response_data = [
            {"population": "AFR", "frequency": 0.5, "n_effective": 10},
            {"population": "EUR", "frequency": 0.3, "n_effective": 8},
        ]
        logging.debug(f"Mock response data: {mock_response_data}")
        
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = mock_response_data

        repeat_id = "123"
        result = GetFreqSTRInfoAPI(repeat_id)

        # Assertions for returned grouped DataFrame
        logging.debug(f"Result from GetFreqSTRInfoAPI: {result}")
        self.assertIsNotNone(result)
        self.assertIn("AFR", result.groups.keys())
        self.assertIn("EUR", result.groups.keys())
        
        logging.info("Completed test: test_GetFreqSTRInfoAPI_success")

    @patch('locus_view.requests.get')
    def test_GetFreqSTRInfoAPI_no_data(self, mock_get):
        """Test GetFreqSTRInfoAPI when no frequency data is returned."""
        logging.info("Starting test: test_GetFreqSTRInfoAPI_no_data")
        
        # Simulate empty API response
        mock_get.return_value.status_code = 200
        mock_get.return_value.json.return_value = []

        repeat_id = "123"
        result = GetFreqSTRInfoAPI(repeat_id)

        # Assertions - should return None for empty response
        logging.debug(f"Result from GetFreqSTRInfoAPI with no data: {result}")
        self.assertIsNone(result)
        
        logging.info("Completed test: test_GetFreqSTRInfoAPI_no_data")

if __name__ == '__main__':
    unittest.main()
