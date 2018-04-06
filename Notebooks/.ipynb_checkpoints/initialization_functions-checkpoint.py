{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "################################\n",
    "# SSH initialization functions #\n",
    "################################\n",
    "\n",
    "\n",
    "def reststate(n_ens=1,path='.')\n",
    "    if n_ens>1:\n",
    "        print('Error: reststate only works for one-member-ensemble')\n",
    "        return\n",
    "    fid = nc.Dataset(path)\n",
    "    lon=np.array(fid.variables[\"nav_lon\"][:])\n",
    "    lat=np.array(fid.variables[\"nav_lat\"][:])\n",
    "    field=np.zeros([1,lon,lat])\n",
    "    field[0,:,:]=np.array(fid.variables[\"degraded_sossheig\"][:,:])\n",
    "    \n",
    "    return field, lon, lat"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
