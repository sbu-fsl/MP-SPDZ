# this tests the pretrained VGG in secure computation

program.options_from_args()

from Compiler import ml

try:
    ml.set_n_threads(int(program.args[2]))
except:
    pass

try:
    batch_size = int(program.args[3])
except:
    batch_size = 1

import torchvision
import torch
import numpy
import requests
import io
import PIL

from torchvision import transforms

name =  'vgg' + program.args[1]
model = getattr(torchvision.models, name)(weights='DEFAULT')

r = requests.get('https://github.com/pytorch/hub/raw/master/images/dog.jpg')
input_image = PIL.Image.open(io.BytesIO(r.content))
input_tensor = transforms._presets.ImageClassification(crop_size=224)(input_image)
input_batch = input_tensor.unsqueeze(0) # create a mini-batch as expected by the model

with torch.no_grad():
    output = int(model(input_batch).argmax())
    print('Model says %d' % output)

input_reshaped = numpy.moveaxis(input_batch.numpy(), 1, -1)

no_input = 'noinput' in program.args

if no_input:
    secret_input = sfix.Tensor([batch_size] + list(input_reshaped.shape[1:]))
else:
    assert batch_size == 1
    secret_input = sfix.input_tensor_via(0, input_reshaped)

layers = ml.layers_from_torch(model, secret_input.shape, batch_size,
                              input_via=None if no_input else 0)

optimizer = ml.Optimizer(layers, time_layers=True)

start_timer(1)
print_ln('Secure computation says %s',
         optimizer.eval(secret_input, top=True)[0].reveal())
stop_timer(1)
