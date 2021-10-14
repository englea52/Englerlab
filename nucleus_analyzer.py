#!/usr/local/bin/python3
"""This program analyzes spheroids!

Please adhere to the Google Python Style Guide: https://google.github.io/styleguide/pyguide.html

Check for style and syntax errors using: "pylint spheroid_analyzer.py"
"""

import argparse
import csv
import concurrent.futures
import math
import os
from os import path
import traceback

import cv2
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy
from PIL import Image
from itertools import zip_longest


DIST = [50, 100, 150, 200, 250, 300, 350, 400]
INT_0 = [1.00, 0.98, 0.95, 0.94, 0.05, 0.05, 0.01, 0.01]
INT_48 = [1.00, 0.98, 0.70, 0.85, 0.95, 0.30, 0.10, 0.01]
MICRONS_PER_PIXEL = 0.09895


def __parse_args():
    """Parses command line arguments
    Example: parser.add_argument("--t2-images", help="t2 image directory",
        required=True, type=str)

    --t2-images: variable name that the input will be assigned to
    help: what input should be provided for the variable name
    required: if the input is required for the program to proceed, required will be True
    type: the input argument type
    """
    parser = argparse.ArgumentParser(description="Generates graphs for spheroid experiments.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--t2-images", help="image directory", required=True, type=str)
    parser.add_argument("--t2-com", help="center of mass file", required=True, type=str)
    parser.add_argument("-o", "--output-dir", default=".", help="Path to output directory.",
                        required=False, type=str)
    #default="." will save the output in the current directory if no output directory is
    #provided
    args = parser.parse_args()
    return args


def __file_listmaker(directory):
    """Makes separate lists of red images and green images in a directory"""
    files = []
    #os.listdir(directory) gets list of files and directories in "directory"
    for filename in os.listdir(directory):
        filepath = path.join(directory, filename)
        if filename[:2] != "._":
            #Ignoring metadata files if source folder is on external hard drive
            print(filename)
            #change for file type and base name
            if "wk" in filename and filename.endswith(".tif"):
                files.append(filepath)
            # if "TXRD" in filename and filename.endswith(".tif"):
            #     red_files.append(filepath)
            # if "FITC" in filename and filename.endswith(".tif"):
            #     green_files.append(filepath)
        else:
            print("Ignoring file: {}".format(filepath))
    return sorted(files)


def __get_bucket(pixel, com, bucket_interval_pixels):
    """
    Returns an index indicating which bucket the pixel falls in. pixel is a pair
    of x-y coordinates. com is a pair of x-y coordinates. bucket_interval_pixels
    is a measure of the difference in length of a bucket's inner and outer
    radii.
    """

    #distance from center of mass (com)
    #numpy.linalg.norm calculates distances between two points
    distance_from_com = numpy.linalg.norm([x - y for x, y in zip(pixel, com)])
    #figuring out which bucket the pixel falls in
    return math.floor(distance_from_com/float(bucket_interval_pixels))


def __get_wedge(pixel, com, wedge_interval_radians):
    """
    Returns an index indicating which wedge the pixel falls in. pixel is a pair
    of x-y coordinates. wedge_interval_radians is the angle that a wedge spans.
    """
    #calculate angle between 0 degree line and the pixel. If the angle is within
    #a given range, then the pixel is in that wedge
    #Formula for calculating angle between pixel and 0 radian point
    angle_rad = math.atan2(pixel[1] - com[1], pixel[0] - com[0])
    if angle_rad < 0:
        angle_rad += 2 * math.pi
    return math.floor(angle_rad/float(wedge_interval_radians))


def __analyze_image(image, com, bucket_interval_pixels, wedge_interval_radians):
    """Analyzes average pixel intensity as a function of center of mass and
    pixel intensity in wedges
    """
    print("Analyzing image: {}".format(image))
    image = Image.open(image)
    pixels = image.load()
    #x pixels = 1392; y pixels = 1040
    width, height = image.size
    min_pixel_value = numpy.amin(image)
    #max_pixel_value = numpy.amax(image)
    if com[0] < width/2.0:
        if com[1] < height/2.0:
            furthest_corner = (width - 1, height - 1)
        else:
            furthest_corner = (width - 1, 0)
    else:
        if com[1] < height/2.0:
            furthest_corner = (0, height - 1)
        else:
            furthest_corner = (0, 0)
    max_distance_from_com = numpy.linalg.norm([x - y for x, y in zip(furthest_corner, com)])
    # Number of concentric circles
    num_buckets = int(math.ceil(max_distance_from_com/float(bucket_interval_pixels)))
    # Number of wedges image will be divided into
    num_wedges = int(math.ceil((2*math.pi)/float(wedge_interval_radians)))

    # List of wedges, where each wedge is a list of buckets. Each bucket is a tuple where the first
    # element is an intensity sum and the second element is a pixel counter.
    wedges_data = [[[0, 0] for _ in range(num_buckets)] for _ in range(num_wedges)]
    target_total_pixels = width * height
    for x_coor in range(width):
        for y_coor in range(height):
            pixel = (x_coor, y_coor)
            bucket_idx = __get_bucket(pixel, com, bucket_interval_pixels)
            wedge_idx = __get_wedge(pixel, com, wedge_interval_radians)

            wedges_data[wedge_idx][bucket_idx][0] += pixels[pixel]
            wedges_data[wedge_idx][bucket_idx][1] += 1

    total_pixels = 0
    for wedge in wedges_data:
        for bucket in wedge:
            total_pixels += bucket[1]
    assert total_pixels == target_total_pixels, "Incorrect number of pixels found!"

    buckets_data = [(sum([x for x in zip(*wedges)][0]), sum([s for s in zip(*wedges)][1]))
                    for wedges in zip(*wedges_data)]
    return buckets_data, wedges_data, min_pixel_value


def __coms_csv_reader(coms_file):
    """Opens coms csv files as lists"""
    with open(coms_file, encoding='utf-8-sig') as csv_file:
        return [(int(row[0]), int(row[1])) for row in csv.reader(csv_file)]


def __nuclear_edge(y_values):
    radius = 0
    for i, y_value in enumerate(y_values):
        if y_value <= 0.2:
            radius = i
            break
    #print("Area:{}".format(math.pi * (radius_t2 ** 2)))
    return radius


def __new_spheroid_radius(image, center, new_radius, output_dir):
    img = cv2.imread(image)
    radius_in_pixels = math.ceil(new_radius / MICRONS_PER_PIXEL)
    cv2.circle(img, center, radius_in_pixels, color=(255, 255, 0), thickness=4)
    filename_prefix = path.basename(image).split(".")[0]
    cv2.imwrite(path.join(output_dir, "{}_new_boundary.png".format(filename_prefix)), img)


def __max_intensity_xdistance(xvalues, yvalues, edge_of_wedge):
    for i, y in enumerate(yvalues):
        if y == 1:
            return xvalues[i]/edge_of_wedge
        else:
            continue


def __wedge_edge(zipped_list):
    for x, dydx in zipped_list:
        if dydx >= 0:
            return x
        else:
            continue


def __bucket_analyzer(t2_file, t2_com, bucket_interval_pixels, wedge_interval_radians,
                      output_dir):
    """analyze pixel intensity for image at t2, some time after initial implanting"""

    t2_analysis, t2_wedge_analysis, _ = __analyze_image(t2_file, t2_com, bucket_interval_pixels,
                                                        wedge_interval_radians)
    x_values = []
    y_values_t2 = []
    bucket_start_pixel = 0
    for (t2_intensity, t2_count) in t2_analysis:
        x_values.append(bucket_start_pixel * MICRONS_PER_PIXEL)
        y_values_t2.append(t2_intensity / float(t2_count))
        bucket_start_pixel += bucket_interval_pixels
    max_y_t2, min_y_t2 = max(y_values_t2), min(y_values_t2)
    y_values_t2 = [(val_t2 - min_y_t2) / (max_y_t2 - min_y_t2) for val_t2 in y_values_t2]
    filename_prefix = path.basename(t2_file).split(".")[0]
    centroid, max_I, ratio = __monochrome_wedge_analyzer(t2_wedge_analysis, bucket_interval_pixels,
                                                         filename_prefix,
                                                         output_dir, wedge_interval_radians)
    #__new_spheroid_radius(t2_file, t2_com, x_values[__nuclear_edge(y_values_t2)], output_dir)
    __wedges_over_image(t2_file, t2_com, filepath_prefix=path.join(output_dir, filename_prefix),
                        wedge_interval_radians=wedge_interval_radians)
    #print("Boundary: {}".format(x_bound))
    return centroid, max_I, ratio, filename_prefix


def __monochrome_wedge_analyzer(wedge_analysis, bucket_interval_pixels, filename_prefix,
                                output_dir, wedge_interval_radians):
    wedge_index = 1
    normalized_wedges = []
    #wedge_derivatives_data = []
    #derivatives_output = []
    centroid_intensity = [('Wedge Index', 'Normalized Intensity')]
    max_intensity_distance = [('Wedge Index', 'Max Intensity Distance from Centroid')]
    intensity_ratios = [('Wedge Index', 'Edge Distance', 'Ratio of Edge to Centroid Intensity')]
    for wedge in wedge_analysis:
        x_values = []
        y_values = []
        bucket_start_pixel = 0
        for (intensity, count)in wedge:
            if count == 0:
                continue
            x_values.append(bucket_start_pixel * MICRONS_PER_PIXEL)
            y_values.append(intensity / float(count))
            bucket_start_pixel += bucket_interval_pixels
        max_y, min_y = max(y_values), min(y_values)
        y_values = [(val - min_y) / (max_y - min_y) for val in y_values]
        nuclear_radius_index = __nuclear_edge(y_values)
        x_values = [x for x in x_values[:nuclear_radius_index]]
        y_values = [y for y in y_values[:nuclear_radius_index]]

        centroid_intensity.append((wedge_index, y_values[0]))
        der_yvalues = numpy.diff(y_values)/numpy.diff(x_values)
        edge_of_wedge = __wedge_edge(zip(reversed(x_values[1:]), reversed(der_yvalues)))
        max_intensity_distance.append((wedge_index, __max_intensity_xdistance(x_values, y_values,
                                                                              edge_of_wedge)))
        edge_index = x_values.index(edge_of_wedge)
        intensity_ratios.append((wedge_index, edge_of_wedge, y_values[edge_index]/y_values[0]))
        normalized_wedges.append([[x, y] for x, y in zip(x_values, y_values)])
        #wedge_derivatives_data.append([[x, y] for x, y in zip(x_values[1:], der_yvalues)])
        #derivatives_output.append([d for d in der_yvalues])

        wedge_index += 1

    __plot_wedges(normalized_wedges, filepath_prefix=path.join(output_dir, filename_prefix),
                  wedge_interval_radians=wedge_interval_radians)
    __plot_derivatives(wedge_derivatives_data, filepath_prefix=path.join(output_dir,
                                                                         filename_prefix),
                       wedge_interval_radians=wedge_interval_radians)

    #derivatives_transpose = [list(row) for row in list(zip_longest(*derivatives_output,
    #                         fillvalue=0))]

    # with open(os.path.join(output_dir, "{}_Wedges Derivatives.csv".format(filename_prefix)),
    #           "w") as wedge_derivatives:
    #     writer = csv.writer(wedge_derivatives, lineterminator="\n")
    #     writer.writerows(derivatives_transpose)

    with open(os.path.join(output_dir, "{}_Centroid Intensities.csv".format(filename_prefix)),
          "w") as centroid_intensities:
        writer = csv.writer(centroid_intensities, lineterminator="\n")
        writer.writerows(centroid_intensity)

    with open(os.path.join(output_dir, "{}_Max Intensity Distance.csv".format(filename_prefix)),
          "w") as max_intensities:
        writer = csv.writer(max_intensities, lineterminator="\n")
        writer.writerows(max_intensity_distance)

    with open(os.path.join(output_dir,
                           "{}_Ratio of Edge to Centroid Intensity.csv".format(filename_prefix)),
                           "w") as ratio:
        writer = csv.writer(ratio, lineterminator="\n")
        writer.writerows(intensity_ratios)

    centroid_intensity_list = [j for i, j in centroid_intensity[1:]]
    centroid_intensity_average = sum(centroid_intensity_list)/len(centroid_intensity_list)

    max_intensity_list = [j for i, j in max_intensity_distance[1:]]
    max_intensity_average = sum(max_intensity_list)/len(max_intensity_list)

    intensity_ratios_list = [k for i, j, k in intensity_ratios[1:]]
    intensity_ratio_average = sum(intensity_ratios_list)/len(intensity_ratios_list)

    return centroid_intensity_average, max_intensity_average, intensity_ratio_average

def  __analyze_files(t2_images, t2_com, output_dir):
    """Goes through the files in directories, finds pictures, processes them
    """
    #Separate red and green files in the given directory
    t2_red = __file_listmaker(t2_images)
    #open the csv file as a list of x,y coordinates
    t2_coms = __coms_csv_reader(t2_com)
    #width of each 'bucket', or concentric circle
    bucket_interval_pixels = 1
    #dividing number is the number of wedges you want to break it into
    wedge_interval_radians = math.pi / 8
    images_data = [('File Name', 'Centroid Intensity', 'Max Intensity Distance',
                   'Edge to Centroid Intensity Ratio')]

    print('{}: Number of images; {}: Number of com values'.format(len(t2_red), len(t2_coms)))
    if len(t2_red) != len(t2_coms):
        raise Exception("Must supply the same number of image files and center-of-mass values.")

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        if t2_red != []:
            print("Performing bucket analysis...")
            analyze_r_futures = [executor.submit(__bucket_analyzer, red_2,
                                                 t2_com, bucket_interval_pixels,
                                                 wedge_interval_radians, output_dir)
                                 for red_2, t2_com in zip(t2_red, t2_coms)]
            for future in concurrent.futures.as_completed(analyze_r_futures):
                try:
                    centroid, max_I, ratio, filename_prefix = \
                    future.result()
                    images_data.append((filename_prefix, centroid, max_I, ratio))

                except Exception as exc:
                    print('{} generated an exception: {}'.format(future, exc))
                    print(traceback.format_exc())

    with open(os.path.join(output_dir, "Results.csv"), "w") as compiled_data:
        writer = csv.writer(compiled_data, lineterminator="\n")
        writer.writerows(images_data)


def __plot_intensity(dist, int_2, filepath_prefix):
    """Generates plots of average intensity as a function of center of mass of spheroid

    Overlays t = 0hrs and t = 48hrs on same graph.

    dist: vectors of distance from center of mass.
    int: vectors of average intensity at distance x from center of mass.
    """
    _, axes = plt.subplots()
    #plot line
    axes.plot(dist, int_2, linestyle="--", color="purple", linewidth=2, label="t = 24 hours")
    #boundary line
    axes.axhline(0.025, color='k', linestyle='--')
    #figure labels
    #axes.legend(loc="upper right")
    axes.set_title("Average Fluorescence Intensity vs. Distance")
    axes.set_xlabel("Distance from center-of-mass ("u"\u03BCm)")
    axes.set_ylabel("Average Intensity\n(Normalized by Max Pixel Intensity)")
    plt.savefig("{}_OverallIntensityPlot.pdf".format(filepath_prefix))


def __plot_wedges(wedges_plot_data, filepath_prefix, wedge_interval_radians):
    """Generates one subplot per wedge comparing red and green intensity at final time"""
    plot_index = 1
    rows = int(math.ceil(len(wedges_plot_data) / 2))
    columns = 2
    fig = plt.figure(figsize=[8, 11])
    for wedge in wedges_plot_data:
        xvalues = []
        yvalues = []
        for x_val, y_val in wedge:
            xvalues.append(x_val)
            yvalues.append(y_val)
        # der_yvalues = numpy.diff(yvalues)/numpy.diff(xvalues)
        ax1 = plt.subplot(rows, columns, plot_index)
        ax1.plot(xvalues, yvalues, color="purple", linewidth=2, label="I(x)")
        ax1.legend(loc="upper right", fontsize=6)
        ax1.set_title("Wedge {}: {} to {} radians".format(plot_index, round((plot_index - 1) * \
                      wedge_interval_radians, 2), round(plot_index * wedge_interval_radians, 2),
                                                          fontsize=7))
        plot_index += 1
    fig.suptitle("Average Fluorescence Intensity vs. Distance")
    fig.text(0.5, 0.005, "Distance from center-of-mass ("u"\u03BCm)", ha='center')
    fig.text(0.0001, 0.5, "Average Intensity\n(Normalized by Max Pixel Intensity", va='center',
             rotation='vertical')
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])
    plt.savefig("{}_WedgesSubplots.pdf".format(filepath_prefix), bbox_inches="tight")

def __plot_derivatives(wedge_derivatives_data, filepath_prefix, wedge_interval_radians):
    plot_index = 1
    rows = int(math.ceil(len(wedge_derivatives_data) / 2))
    columns = 2
    fig = plt.figure(figsize=[8, 11])
    for wedge in wedge_derivatives_data:
        xvalues = []
        yvalues = []
        for x_val, y_val in wedge:
            xvalues.append(x_val)
            yvalues.append(y_val)
        # der_yvalues = numpy.diff(yvalues)/numpy.diff(xvalues)
        ax1 = plt.subplot(rows, columns, plot_index)
        ax1.plot(xvalues, yvalues, color="purple", linestyle="--", linewidth=2, label="dI/dx")
        ax1.legend(loc="upper right", fontsize=6)
        ax1.set_title("Wedge {}: {} to {} radians".format(plot_index, round((plot_index - 1) * \
                      wedge_interval_radians, 2), round(plot_index * wedge_interval_radians, 2),
                                                          fontsize=7))
        plot_index += 1
    fig.suptitle("Derivative of Fluorescence Intensity vs. Distance")
    fig.text(0.5, 0.005, "Distance from center-of-mass ("u"\u03BCm)", ha='center')
    fig.text(0.0001, 0.5, "Derivative of Average Intensity\n(Normalized by Max Pixel Intensity",
             va='center', rotation='vertical')
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])
    plt.savefig("{}_DerivativesWedgeSubplots.pdf".format(filepath_prefix), bbox_inches="tight")

def __wedges_over_image(image, com, filepath_prefix, wedge_interval_radians):
    """Generates plots of average intensity as a function of center of mass of spheroid

    Overlays t = 0hrs and t = 48hrs on same graph.

    dist: vectors of distance from center of mass.
    int: vectors of average intensity at distance x from center of mass.
    """
    num_wedges = int(math.ceil((2*math.pi)/float(wedge_interval_radians)))
    wedge_index = 0
    wedge_radians = []
    img = plt.imread(image)
    _, axes = plt.subplots()
    #plot line
    axes.imshow(img)
    #axes.plot(dist, int_2, linestyle="--", color="purple", linewidth=2, label="t = 24 hours")
    #boundary line
    for i in range(num_wedges):
        xvalues = [com[0], com[0] + 50 * math.cos(wedge_index * wedge_interval_radians)]
        yvalues = [com[1], com[1] + 50 * math.sin(wedge_index * wedge_interval_radians)]
        axes.plot(xvalues, yvalues, color=numpy.random.rand(3,),
                  label="Wedge {}".format(wedge_index + 1))
        wedge_index += 1
    axes.axhline(com[1], color='purple', linestyle='--')
    axes.axvline(com[0], color='purple', linestyle='--')
    #figure labels
    axes.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    axes.set_title("Nucleus with wedges overlay")
    axes.set_xlabel("x distance (pixels)")
    axes.set_ylabel("y distance (pixels)")
    plt.savefig("{}_SectionedNucleus.pdf".format(filepath_prefix))

def __main():
    """First function called when program is run, instead of just going from top to bottom

    Directs the code flow

    Must include if __name__ == "__main__"
                     __main()
    """
    args = __parse_args()
    __analyze_files(args.t2_images, args.t2_com, args.output_dir)


if __name__ == "__main__":
    __main()
