import gzip
import time

def log_message(msg):
    """Print timestamped log message"""
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

def read_gzipped_file(filepath):
    """Read gzipped or regular file"""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def parse_chromosome(chr_str):
    """Parse chromosome string to integer"""
    if chr_str.startswith('chr'):
        chr_str = chr_str[3:]
    try:
        return int(chr_str)
    except ValueError:
        return None  # Skip non-numeric chromosomes
    
