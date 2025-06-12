import os
import sys
import shutil

def remove_all_contents(directory):
    """
    Removes all files and subdirectories in the specified directory.
    """
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    try:
        for item in os.listdir(directory):
            item_path = os.path.join(directory, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
                print(f"Removed file: {item_path}")
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
                print(f"Removed directory: {item_path}")
            else:
                print(f"Skipped (unknown type): {item_path}")
    except Exception as e:
        print(f"Error while removing contents: {e}")
        sys.exit(1)

# Example usage
# if __name__ == "__main__":
#     target_directory = input("Enter the directory to clean: ").strip()
#     remove_all_contents(target_directory)

remove_all_contents('./calc_grid_all')
remove_all_contents('./kin_grid')