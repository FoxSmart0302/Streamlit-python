# Importing necessary libraries
import os
import io
import fitz
import shutil
import docx2txt

from PIL import Image
from pdf2docx import image as pdf_image
from PyPDF2 import PdfReader
from pdf2docx import Converter

# Function to get the total number of pages in a PDF file
def get_total_pages(pdf_file):
    with open(pdf_file, 'rb') as pdf:
        pdf_reader = PdfReader(pdf)
        return len(pdf_reader.pages)

# Function to extract images from the last page of a PDF
def extract_images_from_pdf(pdf_file):
    pdf_document = fitz.open(pdf_file)
    total_pages = get_total_pages(pdf_file)
    
    # Load the last page and extract images
    last_page = pdf_document.load_page(total_pages - 1)
    images_extractor = pdf_image.ImagesExtractor.ImagesExtractor(last_page)
    return images_extractor.extract_images()

# Function to save images as TIFF format (skipping the first two images)
def save_images_as_tiff(images):
    for i, image_data in enumerate(images):
        if i == 0 or i == 1:
            continue
        else:
            print(image_data['width'])
            print(image_data['height'])
            # Open the image from bytes and save it as TIFF
            image_bytes = image_data['image']
            image_obj = Image.open(io.BytesIO(image_bytes))
            image_obj.save(f'image_{i}.tif', 'TIFF')

# Function to convert a PDF file to DOCX
def convert_pdf_to_docx(pdf_file, docx_file):
    total_pages = get_total_pages(pdf_file)
    
    # Convert the PDF to DOCX starting from the last page
    converter = Converter(pdf_file)
    converter.convert(docx_file, start=total_pages-1)
    converter.close()

# Function to process DOCX, extract text and remove the DOCX file
def process_docx_to_text_and_images(docx_file, output_folder):
    text = docx2txt.process(docx_file, output_folder)
    os.remove(docx_file)

# Function to remove duplicate images in a folder based on filename numbering
def remove_duplicates_from_folder(folder):
    max_number = -1

    # Find the maximum numbered image file
    for filename in os.listdir(folder):
        if filename.startswith('image') and filename.endswith('.png'):
            try:
                number = int(filename[5:-4])
                if number > max_number:
                    max_number = number
            except ValueError:
                continue

    # Remove all images except the one with the maximum number
    if max_number >= 0:
        for filename in os.listdir(folder):
            if filename.startswith('image') and filename.endswith('.png'):
                try:
                    number = int(filename[5:-4])
                    if number != max_number:
                        file_to_remove = os.path.join(folder, filename)
                        os.remove(file_to_remove)
                except ValueError:
                    continue

# Main function
def main():
    # Input and output file paths
    pdf_file = 'pdfs/48_BOE-A-2023-3087.pdf'
    docx_file = 'docxs/48_BOE-A-2023-3087.docx'
    images_folder = f"images/{pdf_file.split('/')[-1].split('.')[0]}"

    # Extract images from PDF and save as TIFF
    images = extract_images_from_pdf(pdf_file)
    save_images_as_tiff(images)

    # Convert PDF to DOCX and process the DOCX file
    convert_pdf_to_docx(pdf_file, docx_file)
    process_docx_to_text_and_images(docx_file, images_folder)

    # Remove duplicate images from the output folder
    remove_duplicates_from_folder(images_folder)

# Execute the main function if the script is run directly
if __name__ == "__main__":
    main()
